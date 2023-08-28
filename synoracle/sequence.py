"""
A module containing the tools to analyse synthesis Sequences and produce useful data structures for further analysis.

Author: Joe Manning (@jrhmanning, joseph.manning@manchester.ac.uk)
Date: Nov 2022

Classes:
Sequence - A class containing the linear set of actions within a synthesis and methods for analysing them

Exceptions:
None

TODO: improve the prioritise_bom_df and produce_bill_of_mats workflow to reduce number of possible errors
TODO: output ingredients.BillOfMaterials from produce_bill_of_mats
"""
from __future__ import annotations
import json # not needed?
import pandas as pd
import logging
import numpy as np
import re # not needed?
from pint import UndefinedUnitError #not needed?
from pint import UnitRegistry, Unit
from collections import defaultdict
from typing import Dict
from .ingredients import ChemicalList
from .conditions import Conditions
from textdistance import levenshtein

ureg = UnitRegistry()
Q_ = ureg.Quantity

# create logger
logger = logging.getLogger('simple_example.txt')
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)

class SequenceGraph():
    def __init__(self):
        pass


class Sequence():
    """
    A class for handling a synthesis sequence as a table, including methods to eliminate spurious data and sumarise
    overly verbose explanations into more information-dense alternatives.
    As the quesence is primarily a pandas DataFrame, the majority of methods just modify specific aspects of the df.

    Key methods:
    : drop_invalid_lines: removes any empty actions from the sequence
    : assign_step_supertypes: provides a class of each action, reducing ChemicaTagger's 15 action types to 3-4
    : condense_to_supertypes: Summarises data from adjacent actions of a specific supertype to a single action
    : extract_chemicals: pulls out a flattened dataframe of chemicals mentioned in the sequence
    : prioritise_bom_df: associates each chemical mention with at most one quantity (prevents double-counting)
    : produce_bill_of_mats: aggregates a chemical list into a condensed bill of materials (defunct?)
    : parse_temperature_to_kelvin: Uses pint to standardise temperature mentions from strings to floats in Kelvin
    : parse_time_to_minutes: Uses pint to standardise time mentions from strings to floats in minutes
    : extract_times_and_temps: Applies the above two functions to the whole sequence

    Static methods
    : from_json: instantiates from a json object, to help interface with databases

    """
    type_subs = {  # key = ChemicalTagger action type; value = My abbreviation of it
        'add': 'add',  # add new chemicals
        'apparatusaction': 'react',  # signifies a reaction happening
        'concentrate': 'condition',  # signifies a non-reactive step
        'cool': 'condition',
        'degass': 'remove',  # signifies a workup/postprocessing step
        'dissolve': 'add',
        'dry': 'remove',
        'extract': 'remove',
        'filter': 'remove',
        'heat': 'condition',
        'partition': 'remove',
        'precipitate': 'remove',
        'purify': 'remove',
        'quench': 'remove',
        'recover': 'remove',
        'remove': 'remove',
        'stir': 'add',
        'synthesize': 'react',
        'wait': 'react',
        'wash': 'remove',
        'yield': 'remove'
    }

    def __init__(self, raw_sequence: pd.DataFrame):
        """ Instantiates the object from a pandas DataFrame"""
        self.raw_synthesis = raw_sequence.rename({'Step name': 'name'}, axis=1)
        # TODO: unit test for handling different data types as raw_sequence
        self.drop_invalid_lines()

    def drop_invalid_lines(self):
        """
        Deletes all actions with no chemical, time, or temperature quantities included
        :return: None
        """
        # TODO: unit tests to confirm that lines are being dropped correctly
        def test_empty(df_cell):
            """
            Tests if a quantity is present in a cell of the raw sequence.
            Valid input data types are floats, strings, lists, dicts, sets, or 
            strings that turn into lists etc. when called with eval()). 
            Returns True if the input data is not np.nan or an empty iterable
            """
            # TODO: unit tests to confirm different cell data are handled correctly 

            try: # turn df-safe python data structures into their structures
                working = eval(df_cell)
            except SyntaxError: # means it's a string(?)
                if len(df_cell)>1:
                    logging.debug(f'String found during drop_invalid_lines: {df_cell}')
                    return True
            except TypeError: # means it's either an iterable or np.nan
                try:
                    iter(df_cell)
                except TypeError: # means it's a float
                    return df_cell is not np.nan
                else: # it's an iterable!
                    return len(df_cell)>0
            else:
                return len(working)>0
            raise TypeError("Couldn't parse dataframe cell information: {0}".format(df_cell))

        new_df = pd.DataFrame(columns=self.raw_synthesis.columns)
        for c,row in self.raw_synthesis.iterrows():
            cont_chems = test_empty(row['new_chemicals'])
            cont_times = test_empty(row['time'])
            cont_temps = test_empty(row['temp'])
            if any([cont_chems, cont_times, cont_temps]):
                new_df.loc[c, :] = row

        self.clean_synthesis = new_df.reset_index(drop=True)
        logging.debug(self.clean_synthesis)

    def assign_step_supertypes(self):
        """
        Categorises sequence data by action "supertype" by a similar method to the ULSA paper by Kononova
        Directly maps the 15 action types from ChemicalTagger to 3-4 types, and drops "Start" statements
        :return:
        """
        # TODO: Allow multiple start statements at the beginning of the sequence
        # TODO: unit tests to check if mapping is correctly working
        # TODO: unit tests to check handling of invalid data mappings

        try:
            self.clean_synthesis
        except AttributeError:
            self.drop_invalid_lines()

        self.clean_synthesis['Step supertype'] = self.clean_synthesis['name'].apply(
            lambda x: str(x).lower()).map(self.type_subs, na_action='ignore')

        self.clean_synthesis['Step supertype'] = self.clean_synthesis['Step supertype'].fillna('other')

        
        try:                             
            first_add = np.where(self.clean_synthesis['Step supertype']=='add')[0][0]
        except IndexError:
            pass
        else:
            self.clean_synthesis.loc[
                (self.clean_synthesis['step number'] < first_add) & (self.clean_synthesis['Step supertype'].apply(
            lambda x: str(x).lower()).isin(['synthesize', 'react', 'yield'])),
                'Step supertype'
            ] = 'start'

    def _ranges(self,nums):
        nums = sorted(set(nums))
        gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s + 1 < e]
        edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
        return list(zip(edges, edges))

    def find_unit_operations(self, ignore_conditions = True):
        """
        Locates groups of sequential steps within a synthesis sequence, which we term as "unit operations"
        To let us ignore certain types of step, we analyse based on the dataframe index, but return based on the "step number" column.
        From this we can group together context-dependent steps like heating and cooling with their neighbouring context-independent steps like synthesis or drying
        As a caveat, we will currently lose some data due to context-dependent steps being between groups of context-independent steps (rather than within a group of context-independenx steps.)  
        :return: a dict of "unit op type": [unit op positions]
        """
        supertypes = list(set(self.type_subs.values()))

        try: 
            self.clean_synthesis['Step supertype']
        except KeyError:
            self.assign_step_supertypes()

        if ignore_conditions:
            working_sequence = self.clean_synthesis.drop(
                self.clean_synthesis[self.clean_synthesis['Step supertype']=='condition'].index
                ).reset_index()
        else:
            working_sequence = self.clean_synthesis

        groups = defaultdict(list)
        for supertype in supertypes:
            data = list(working_sequence[working_sequence['Step supertype'] == supertype].index.astype(int))

            #groups[supertype].extend(self._ranges(data))
            for x in self._ranges(data):
                groups[supertype].append(
                    (
                        working_sequence.loc[x[0],'step number'],
                        working_sequence.loc[x[1],'step number'],
                        )
                )

        return groups

    def condense_to_supertypes(self, ignore_conditions = True):
        """
        Reduces a synthesis sequence to a sequence of unit operations (rather than actions)
        WARNING: not tested!
        :return: None
        """
        try:
            self.clean_synthesis
        except AttributeError:
            self.drop_invalid_lines()

        try:
            self.clean_synthesis['Step supertype']
        except KeyError:
            self.assign_step_supertypes()

        self.clean_synthesis.index = self.clean_synthesis.index.astype(int)


        supertypes = list(set(self.type_subs.values()))

        groups = self.find_unit_operations()

        ordered_groups = sorted([x for y in list(groups.values()) for x in y], key=lambda x: list(x)[0]) 
          #maps ordered_groups tot he dataframe_Sindex. WARNING: assumes that step numebr monotonically increases!
        condensed_sequence = pd.DataFrame()
        for group in ordered_groups:
            ordered_indices = [np.where(self.clean_synthesis['step number']==x)[0][0] for x in group]
            # print(ordered_indices[0] in self.clean_synthesis.index) # this needs a unit test!
            # TODO: check this is collecting the correct step indices
            logging.debug(group)
            try:
                to_condense = self.clean_synthesis.loc[ordered_indices[0]:ordered_indices[1],
                    ['name', 'Step supertype', 'new_chemicals', 'temp', 'time']
                ].groupby('Step supertype').aggregate(sum)
            except KeyError:
                raise
                # to_condense = self.clean_synthesis.loc[str(group[0]):,
                #     ['Step supertype', 'new_chemicals', 'temp', 'time']
                # ].groupby('Step supertype').aggregate(sum)
            except:
                raise
                # to_condense = self.clean_synthesis.loc[str(group[0]):,:]

            if to_condense.empty:
                assert group[0] in self.clean_synthesis.index,"Step indices don't match up with the dataframe indices - chack if you need to run df.reset_index() first!"
                # TODO: modify the code so
                print(ordered_groups)
                print(group)
                print(to_condense)
                try:
                    print(self.clean_synthesis.loc[group[0],:])
                    print(self.clean_synthesis.loc[group[1],:])
                except KeyError:
                    print(self.clean_synthesis.index, self.clean_synthesis.index.dtype)
                    raise
                raise IndexError('Empty dataframe found when steps were expected')
                # to_condense = self.clean_synthesis.loc[str(group[0]):,:]

            try:
                to_condense.loc[:, 'Condensed steps'] = group[1] - group[0] + 1
            except ValueError:
                logging.debug(to_condense)
                raise
            condensed_sequence = condensed_sequence.append(to_condense)

        self.condensed_sequence = condensed_sequence

    def extract_chemicals(self, partial_sequence = None):
        """
        Extracts all chemicals from a sequence dataframe into a Chemical List dataframe (or object) for further analysis
        TODO: refactor ChemicalList object to be what this points to
        TODO: refactor handling of chemicals from this object into ChemicalList
        TODO: create classmethods so that I can instantiate a ChemicalList from a networkx node attribute
        TODO: test functionality to handle a partial chemical list
        :return: None
        """

        if partial_sequence is None:
            try:
                self.clean_synthesis
            except AttributeError:
                self.drop_invalid_lines()
            sequence_to_be_processed = self.clean_synthesis
        else:
            sequence_to_be_processed = partial_sequence

        assert 'new_chemicals' in sequence_to_be_processed, 'Missing chemical informaiton in dataframe! (column name should be "new_chemicals")'

        self.xml_cems = []
        for _, row in sequence_to_be_processed.iterrows():
            # catch poorly-imported jsons in the sequence
            if isinstance(row['new_chemicals'], str):
                working = eval(row['new_chemicals'])
            else:
                working = row['new_chemicals']
            try:
                self.xml_cems.extend(working)
            except (TypeError, KeyError):
                continue

        output_df = pd.DataFrame(
            columns=['name', 'mass', 'other_amount', 'volume', 'percent', 'concentration', 'aliases']
        )

        output_df = output_df.append(self.xml_cems)

        self.xml_cem_ids = list(set([x['name'] for x in self.xml_cems]))
        self.chemical_list = ChemicalList(
            pd.DataFrame(output_df)
            ) # new version using ChemicalList, for testing

    def extract_conditions(self, partial_sequence = None):
        if partial_sequence is None:
            try:
                self.clean_synthesis
            except AttributeError:
                self.drop_invalid_lines()
            sequence_to_be_processed = self.clean_synthesis
        else:
            sequence_to_be_processed = partial_sequence

        assert 'time' in sequence_to_be_processed, 'Missing "time" column in sequence dataframe!'
        assert 'temp' in sequence_to_be_processed, 'Missing "temp" column in sequence dataframe!'

        self.conditions = Conditions(sequence_to_be_processed)

    def calculate_similarity(self, other: Sequence) -> float:
        """
        Start with Levenshtein distance because I don't know any better.
        """
        return levenshtein.normalized_similarity(
            self.condensed_sequence['Step supertype'].tolist(),
            other.condensed_sequence['Step supertype'].tolist()
         )

    @staticmethod
    def from_json(json_filename):
        working = pd.read_json(json_filename)
        # paper_index = '.'.join([working.loc['0','paper'], working.loc['0','paragraph']])
        # result = working['name'].to_list()
        return Sequence(working)

"""
A module containing the tools to analyse synthesis Sequences and produce useful data structures for further analysis.

Author: Joe Manning (@jrhmanning, joseph.manning@manchester.ac.uk)
Date: Nov 2022

Classes:
Sequence - A class containing the linear set of actions within a synthesis and methods for analysing them
SequenceGraph - empty for now (delete?)

Exceptions:
None

TODO: improve the prioritise_bom_df and produce_bill_of_mats workflow to reduce number of possible errors
TODO: output ingredients.BillOfMaterials from produce_bill_of_mats
"""
from __future__ import annotations
import json
import pandas as pd
import logging
import numpy as np
import re
from pint import UndefinedUnitError
from pint import UnitRegistry, Unit
from collections import defaultdict
from typing import Dict
from .ingredients import ChemicalList, BillOfMaterials
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
        self.drop_invalid_lines()

    def drop_invalid_lines(self):
        """
        Deletes all actions with no chemical, time, or temperature quantities included
        :return: None
        """
        self.clean_synthesis = self.raw_synthesis.dropna(
            how='all',
            subset=['new_chemicals', 'time', 'temp']
        )

        logging.debug(self.clean_synthesis)

    def assign_step_supertypes(self):
        """
        Categorises sequence data by action "supertype" by a similar method to the ULSA paper by Kononova
        Directly maps the 15 action types from ChemicalTagger to 3-4 types, and drops "Start" statements
        TODO: Allow multiple start statements at the beginning of the sequence
        :return:
        """
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
                        working_sequence.loc[x[0],'step number'].astype(int), 
                        working_sequence.loc[x[1],'step number'].astype(int)
                        )
                )
            # groups[supertype].extend([
            #     working_sequence.loc[x,'step number'].astype(int) for x in self._ranges(data)
            #     ])
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
        condensed_sequence = pd.DataFrame()
        for group in ordered_groups:

            # TODO: check this is collecting the correct step indices
            logging.debug(group)
            try:
                to_condense = self.clean_synthesis.loc[group[0]:group[1],
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
                print(ordered_groups)
                print(group)
                print(to_condense)
                raise
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

    #region moved to ingredients.ChemicalList
    # def compress_values(self, col):
    #     """Possibly defunct helper funciton to condense info in dataframes to lists"""
    #     return col[col.notna()].to_list() if len(col[col.notna()].to_list()) > 0 else np.nan

    # def prioritise_bom_df(self):
    #     """
    #     Prevents double counting of quantities for each chemical by identifying the most reliable quantity used for each
    #      chemical mention.

    #     Here's how it works:
    #     First define unit priorities. Higher numbers supersede lower. I can't remember why I picked this specific
    #     priority order, but it seems to work. Most importantly moles > anything else > percent; percent is usually
    #     only used for purity At some point it might be interesting to try reordering concentration, volume,
    #     and mass and checking discrepancies.
    #     Then loop through a dataframe of chemical mentions providing 1 priority per row.
    #     :return: None
    #     """

    #     priorities = {'other_amount': 5,
    #                   'concentration': 4,
    #                   'volume': 3,
    #                   'mass': 2,
    #                   'percent': 1
    #                   }

    #     logging.info(priorities)
    #     priority_key = {}

    #     for rownum, row in self.chemical_list.loc[:,
    #                        ['mass', 'concentration', 'other_amount', 'volume', 'percent']
    #                        ].iterrows():

    #         logging.debug(type(rownum))
    #         working_priorities = row[row.notna()].index.values
    #         if len(working_priorities) == 0:
    #             priority_key[rownum] = np.nan
    #             continue

    #         max_priority = max([priorities[x] for x in working_priorities])
    #         priority_key[rownum] = list(priorities.keys())[list(priorities.values()).index(
    #             max_priority)]  # this desperately needs a refactor for neatness, but works(?)

    #     logging.debug(priority_key)
    #     self.chemical_list.loc[:, 'Units used'] = pd.Series(priority_key)

    # def produce_bill_of_mats(self):
    #     """
    #     Aggregates a chemical list by the identified name (for now), aggregating all unit types to  list.
    #     The Units used list then points to the appropriate unit for each mention that's been made.
    #     TODO: decide if it's better to calculate moles first, then aggregate once all the mentions have been standardised
    #     :return: None
    #     """

    #     aggregation_functions = {'name': 'first',
    #                              'mass': self.compress_values,
    #                              'other_amount': self.compress_values,
    #                              'volume': self.compress_values,
    #                              'percent': self.compress_values,
    #                              'concentration': self.compress_values,
    #                              'aliases': self.compress_values,
    #                              #     'pubchem_id': 'first',
    #                              'Units used': self.compress_values
    #                              }

    #     new_ingreds = self.chemical_list.groupby(self.chemical_list['name']).aggregate(
    #         aggregation_functions).reset_index(drop=True)
    #     new_ingreds = new_ingreds.drop('aliases', axis=1)
    #     new_ingreds = new_ingreds.dropna(how='all', axis=1)

    #     #    new_ingreds['pubchem_id'] = new_ingreds['pubchem_id'].astype(int)

    #     self.chemical_list_reduced = new_ingreds

    #endregion

    #region moved to new Conditions class
    # temp_dict = {
    #     '°C': ureg.degC,
    #     'kelvin': ureg.kelvin,
    #     'Kelvin': ureg.kelvin,
    #     'k': ureg.kelvin,
    #     'K': ureg.kelvin
    # }

    # def parse_temperature_to_kelvin(self, df_temp_string: str) -> float:
    #     """
    #     Calculates the numerical temperature from a string of the temperature value + unit
    #     Performs some regex substitutions to recognise common phrases like "room temperature"
    #     :param df_temp_string: a string
    #     :return: The temperature, in kelvin
    #     """
    #     working = re.sub(r'\b(at|of|in|to|for)\b', '', df_temp_string).strip()
    #     try:
    #         output = Q_(float(working.split()[-2]), self.temp_dict[working.split()[-1]])
    #     except (ValueError, KeyError, UndefinedUnitError):
    #         if bool(re.search(r'(R|room|ambient|indoor)\Wtemperature|RT', working)):
    #             output = Q_(25, ureg.parse_units('degC'))
    #         else:
    #             print(f'Failed for: "{working}"')
    #             return np.nan
    #     except:
    #         print(f'Temperature extraction failed for: "{working}", original was "{df_temp_string}"')
    #         return np.nan
    #     # print(output)
    #     return output.to_base_units().magnitude

    # def parse_time_to_minutes(self, df_time_string: str) -> float:
    #     """
    #     Identifies a numerical time from a string of time quantity + unit.
    #     Performs regex to replace most common values for pint-parseable alternatives.
    #     :param df_time_string: a timephrase string
    #     :return: the amount of minutes used as a float
    #     """
    #     working = re.sub(r'\b(at|of|in|to|for|(A|a)fter)\b', '', df_time_string).strip()
    #     working = re.sub(r'~|∼', '', working).strip()
    #     working = re.sub(r'\b(h)\b', 'hour', working).strip()
    #     working = re.sub(r'\b(d|D)\b', 'day', working).strip()

    #     working = re.sub(r'\b(an|a|one)\W(hour|day)\b', r'1 \2', working).strip()
    #     working = re.sub(r'\b(one)\W', r'1 ', working).strip()
    #     working = re.sub(r'\b(two)\W', r'2 ', working).strip()
    #     working = re.sub(r'\b(three)\W', r'3 ', working).strip()
    #     working = re.sub(r'\b(four)\W', r'4 ', working).strip()
    #     working = re.sub(r'\b(five)\W', r'5 ', working).strip()
    #     working = re.sub(r'\b(six)\W', r'6 ', working).strip()
    #     working = re.sub(r'\b(seven)\W', r'7 ', working).strip()
    #     working = re.sub(r'\b(eight)\W', r'8 ', working).strip()
    #     working = re.sub(r'\b(nine)\W', r'9 ', working).strip()

    #     working = re.sub(r'\b(mins)\b', r'minutes ', working).strip()

    #     if bool(re.search('(O|o)vernight|(O|o)ne night', working)):
    #         output = Q_(18, ureg.hour)

    #     else:
    #         try:
    #             output = Q_(float(working.split()[-2]), working.split()[-1])
    #         except (ValueError, KeyError, UndefinedUnitError, IndexError):
    #             print(f'Time extraction failed for: "{working}", original was "{df_time_string}"')
    #             return np.nan
    #     return output.to('minutes').magnitude

    # def extract_times_and_temps(self):
    #     """
    #     Produces a dataframe of quoted times and temperatures for a whole sequence in standardised units.
    #     TODO: add functionality to work for a partial sequence too.
    #     :return: None
    #     """
    #     output = self.clean_synthesis[['time', 'temp']]
    #     for c, x in output.iterrows():
    #         logging.debug(c, x['temp'], x['time'])
    #         if isinstance(x['temp'], str):
    #             x['temp'] = eval(x['temp'])
    #         if isinstance(x['time'], str):
    #             x['time'] = eval(x['time'])
    #         try:
    #             output.loc[c, 'T (K)'] = self.parse_temperature_to_kelvin(x['temp'][0])
    #         except (IndexError, TypeError):
    #             pass
    #         try:
    #             output.loc[c, 'Time (min)'] = self.parse_time_to_minutes(x['time'][0])

    #         except (IndexError, TypeError):
    #             pass

    #     # output['T (K)'] = output['temp'].apply(lambda x: self.parse_temperature_to_kelvin(str(x)))
    #     # output['Time (min)'] = output['time'].apply(lambda x: self.parse_time_to_minutes(str(x)))
    #     self.time_temp = output
    #endregion

    @staticmethod
    def from_json(json_filename):
        working = pd.read_json(json_filename)
        # paper_index = '.'.join([working.loc['0','paper'], working.loc['0','paragraph']])
        # result = working['name'].to_list()
        return Sequence(working)

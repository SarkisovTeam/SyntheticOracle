"""
A module containing the tools to analyse the chemicals used in synthesis, both individually though the ChemicalEntity
 and CheimcalSpecies classes (the former containnig quantity information, the latter cheminformatics identifiers),
 and collectively through the ChemicalList and BillofMaterials object. Here, the ChemicalList object is a raw list of
 chemical mentions from a synthesis, while BillofMaterials is grouped by chemical identities used.

Author: Joe Manning (@jrhmanning, joseph.manning@manchester.ac.uk)
Date: Nov 2022

Classes:
ChemicalList - a list of chemical mentions found in a synthesis, and the raw units associated
BillOfMaterials - a list of chemicals grouped by identity, with summarised units
ChemicalSpecies - identity information of a chemical for cross-referencing databases and cheminformatics
ChemicalEntity(ChemicalSpecies) - Quantity information about a chemical, connected to identity information

Exceptions:
ChemicalNotFoundError - raised when database cross-referencing fails during ChemicalSpecies instantiation
"""
from __future__ import annotations
import pandas as pd
import numpy as np
import pubchempy as pcp
from pubchempy import BadRequestError
from pubchempy import NotFoundError # not needed?
from pint import UnitRegistry, Unit
from urllib.error import URLError

ureg = UnitRegistry()
Q_ = ureg.Quantity
from pint import UndefinedUnitError, DefinitionSyntaxError, DimensionalityError
from pint.util import UnitsContainer
from collections import Counter
import json
import chemicals
import logging
from time import sleep

# create logger
logger = logging.getLogger('simple_example.txt')
logger.setLevel(logging.DEBUG)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
logger.addHandler(ch)



class ChemicalNotFoundError(Exception): pass

class EmptyChemicalList(Exception): pass

class ChemicalList:
    """
    A class describing a raw list of chemicals prior to aggregation by chemical rtype or standardisation of units.
    Contains methods to convert a list of chemicals into a bill of materials by (a) converting string representations
    of units into pint units of moles and aggregating by pubchem ID.

    Key methods:
    : prioritise_bom_df: associates each chemical mention with at most one quantity (prevents double-counting)
    : produce_bill_of_mats: aggregates a chemical list into a condensed bill of materials (defunct?)
    : from_json: instantiates from a dict-like json object, for manual chemical list creation/database imports
    """
    def __init__(self, chemical_list: pd.DataFrame):
        self.chemical_list = chemical_list
        if not len(self.chemical_list)>0:
            logging.warn("It looks like you've tried to make a chemical list with no chemicals in it!")
            raise EmptyChemicalList
        self.prioritise_bom_df()


    def compress_values(self, col) -> list:
        """
        Quick helper function to assist with dataframe.groupby(). Turns a column of lists or nan into a single list.
        Possibly defunct because add does the same thing now?
        Either that or enables mixed list/float/nan --> list
        :param col: a column from a pandas dataframe containing either list or np.nan
        :return: a list or np.nan
        """
        return col[col.notna()].to_list() if len(col[col.notna()].to_list()) > 0 else np.nan

    def prioritise_bom_df(self):
        """
        Prevents double counting of quantities for each chemical by identifying the most reliable quantity used for each
         chemical mention.

        Here's how it works:
        First define unit priorities. Higher numbers supersede lower. I can't remember why I picked this specific
        priority order, but it seems to work. Most importantly moles > anything else > percent; percent is usually
        only used for purity At some point it might be interesting to try reordering concentration, volume,
        and mass and checking discrepancies.
        Then loop through a dataframe of chemical mentions providing 1 priority per row.
        :return: None
        """

        priorities = {'other_amount': 5,
                      'concentration': 4,
                      'volume': 3,
                      'mass': 2,
                      'percent': 1
                      }

        logging.info(priorities)
        priority_key = {}

        # Check all keys are present, otherwise raise an error:
        try:
            self.chemical_list.loc[:,
                           ['mass', 'concentration', 'other_amount', 'volume', 'percent']
                           ]
        except KeyError:
            logging.warn('Missing quantity columns in chemical list!\nMay be indicative of no chemicals being there, or a json/xml parsing error!')
            print('Chemicals in list:',len(self.chemical_list))
            print(self.chemical_list)
            raise


        # Then loop through a dataframe of chemical mentions providing 1 priority per row
        for rownum, row in self.chemical_list.loc[:,
                           ['mass', 'concentration', 'other_amount', 'volume', 'percent']
                           ].iterrows():

            logging.debug(type(rownum))
            working_priorities = row[row.notna()].index.values
            if len(working_priorities) == 0:
                priority_key[rownum] = np.nan
                continue

            max_priority = max([priorities[x] for x in working_priorities])
            priority_key[rownum] = list(priorities.keys())[list(priorities.values()).index(
                max_priority)]  # this desperately needs a refactor for neatness, but works(?)

        logging.debug(priority_key)
        self.chemical_list.loc[:, 'Units used'] = pd.Series(priority_key)

    def produce_bill_of_mats(
            self, 
            identifier_cache_location = None,
            property_cache_location = None):
        """
        Aggregates a chemical list by the identified name (for now), aggregating all unit types to  list.
        The Units used list then points to the appropriate unit for each mention that's been made.
        TODO: decide if it's better to calculate moles first, then aggregate once all the mentions have been standardised
        :return: None
        """
        aggregation_functions = {'name': 'first',
                                 'mass': self.compress_values,
                                 'other_amount': self.compress_values,
                                 'volume': self.compress_values,
                                 'percent': self.compress_values,
                                 'concentration': self.compress_values,
                                 'aliases': self.compress_values,
                                 #     'pubchem_id': 'first',
                                 'Units used': self.compress_values
                                 }

        try: 
            self.chemical_list['Unit_used']
        except KeyError:
            self.prioritise_bom_df()

        new_ingreds = self.chemical_list.groupby(self.chemical_list['name']).aggregate(
            aggregation_functions).reset_index(drop=True)
        new_ingreds = new_ingreds.drop('aliases', axis=1)
        new_ingreds = new_ingreds.dropna(how='all', axis=1)

        #    new_ingreds['pubchem_id'] = new_ingreds['pubchem_id'].astype(int)

        return BillOfMaterials(
            new_ingreds, 
            identifier_cache_location=identifier_cache_location, 
            property_cache_location=property_cache_location)

    @staticmethod
    def from_json(json_pointer):
        json_object = json.loads(json_pointer)

        paper_index = '.'.join([json_object["paper"], json_object["paragraph"]])

        logging.info(f'working on paper {paper_index}!')
        try:  # skip procedure if there's no units identified
            json_object['name']
        except KeyError:
            logging.warning(f'paper {paper_index} has no identified chemicals!')
            return pd.DataFrame()
        bill_of_mats = pd.DataFrame(json_object).fillna(np.nan)
        bill_of_mats["paper_index"] = paper_index
        # bill_of_mats['pubchem_id'] = bill_of_mats['name'].apply(find_pubchem_name)
        #
        # for index, row in bill_of_mats.iterrows():
        #     try:
        #         bill_of_mats.loc[index, 'total_moles'] = calc_moles(row)
        #     except KeyError:
        #         bill_of_mats.loc[index, 'total_moles'] = np.nan

        bill_of_mats = bill_of_mats.set_index("paper_index")

        return bill_of_mats

    # TODO: work out the purpose of this class - is it just to create bills of materials?


class BillOfMaterials:
    """
    A bill of materials for a synthesis (or subset of actions therein). Contains all positively identified chemicals
    and associated amounts therein, standardised as moles.

    """
    def __init__(self, input_dataframe: pd.DataFrame, identifier_cache_location = None, property_cache_location = None):
        self.starting_chemical_list = input_dataframe
        try:
            self.starting_chemical_list['pubchem_id']
        except KeyError:
            self.convert_to_bom(self.starting_chemical_list, identifier_cache_location=identifier_cache_location, property_cache_location=property_cache_location)
        else:
            self.initial_bill_of_materials = self.starting_chemical_list
        

        #TODO: change the aggregate for name to make a list if it's a string, and concat a list if it's a list
        self.bill_of_materials = self.initial_bill_of_materials.reset_index().groupby('pubchem_id').aggregate(
            {'name': lambda x: list(x), 'moles': sum}
        )



    def convert_to_bom(self, initial_dataframe, identifier_cache_location, property_cache_location):
        self.initial_bill_of_materials = pd.DataFrame(columns=['moles', 'pubchem_id'])
        for _, x in initial_dataframe.iterrows():
            try:
                thing = ChemicalEntity(
                    x['name'], 
                    identifier_cache_location=identifier_cache_location, 
                    property_cache_location=property_cache_location)
            except ChemicalNotFoundError:
                continue

            self.initial_bill_of_materials.loc[x['name'], 'pubchem_id'] = thing.pubchem_cid

            thing.calc_moles(x)  # unit convert if needed
            logging.debug(f'Calculated moles for molecule {thing.pubchem_cid}: {thing.moles}')

            try:  # catch if moles are already in the columns list
                self.initial_bill_of_materials.loc[x['name'], 'moles'] = x['moles']
            except KeyError:
                self.initial_bill_of_materials.loc[x['name'], 'moles'] = thing.moles
            except:
                self.initial_bill_of_materials.loc[x['name'], 'moles'] = np.nan
        self.initial_bill_of_materials.index = self.initial_bill_of_materials.index.rename('name')
        logging.debug(
            f'Calculated bill of materials (without aggregation by pubchem number): \n{self.initial_bill_of_materials}')

    def calculate_similar_chemicals(self, other: BillOfMaterials) -> float:
        similar_chemicals = set(self.bill_of_materials.index).intersection(set(other.bill_of_materials.index))
        total_chemicals = set(self.bill_of_materials.index).union(set(other.bill_of_materials.index))
        return 100*(len(similar_chemicals)/len(total_chemicals))

    def calculate_similar_quantities(self, other: BillOfMaterials) -> pd.DataFrame:
        similar_chemicals = set(self.bill_of_materials.index).intersection(set(other.bill_of_materials.index))
        logging.debug(similar_chemicals)
        similar_chemical_df = pd.DataFrame()
        for x in similar_chemicals:
            similar_chemical_df.loc[x, 'percentage quantity difference'] = 100*(self.bill_of_materials.loc[x,'moles'] - other.bill_of_materials.loc[x, 'moles'])/self.bill_of_materials.loc[x,'moles']
            similar_chemical_df.loc[x, 'log quantity difference'] = np.log(other.bill_of_materials.loc[x, 'moles']/self.bill_of_materials.loc[x,'moles'])
        
        return similar_chemical_df


    def calculate_weighted_similarity(self, other: BillOfMaterials) -> float:
        total_chemicals = set(self.bill_of_materials.index).union(set(other.bill_of_materials.index))
        similar_chemical_df = self.calculate_similar_quantities(other)
        sim_index = 0
        for _,x in similar_chemical_df.iterrows():
            log_diff = max((1-0.01*x['log quantity difference']), 0)
            sim_index += log_diff 
        
        return 100*log_diff/total_chemicals



class ChemicalSpecies:
    """
    A class containing all essentials information for uniquely identifying a chemical species and cross-referencing
    it against several important cheminformatics packages. Useful for vital unit transformations e.g. volume -> moles
    by estimating pure component density.

    Key methods:
    : find_pubchem_name: finds the pubchem CID for a compound from the identified alias as a name or formula
    : extract_key_properties_from_pubchem: pulls useful properties to reduce memory load/database calls
    : find_chedl: finds the molecule's CAS for looking up the ChEDL database
    : density_from_chedl: if possible, determines the pure component density using the COSTALD algorithm
    """
    def __init__(
        self, 
        alias=None, 
        identifier_cache_location = None,
        property_cache_location = None
        ):

        self.identifier_cache_location = identifier_cache_location
        self.property_cache_location = property_cache_location
        # assert isinstance(alias, str) or isinstance(alias, int), f'Chemical name is invalid, {alias} ({type(alias)})'

        if not alias:
            pass

        if isinstance(alias, int):
            self.pubchem_cid = alias
        elif isinstance(alias, str):
            try:
                int(alias)
                self.pubchem_cid = int(alias)
            except ValueError:
                self.find_pubchem_name_memoized(alias)
        else: 
            raise ChemicalNotFoundError

        if self.pubchem_cid:
            assert isinstance(self.pubchem_cid, int), f'pubchem id is invalid! {self.pubchem_cid}, {alias}'
            self.extract_key_properties_from_pubchem_memoized()
            self.find_chedl()
        else: 
            raise ChemicalNotFoundError

        try:
            self.CAS
            self.density_from_chedl()
        except (AttributeError, TypeError):
            pass

    def __str__(self):
        output = ''
        try:
            output += f'A chemical called {self.name}\n\n'
            output += 'Attributes:\n--------------------------------\n'
            output += f'Molecular formula: {self.mol_formula} ({self.mol_wt} g/mol)\n'
            output += f'Smiles: {self.smiles}\n'
            output += f'CAS ID: {self.CAS}\n'
            # output += f'COSTALD liquid density: {self.density}'
            output += '\n--------------------------------\n'
        except:
            pass
        return output

    def find_pubchem_name(self, alias: str):
        """
        Database lookup function for the pubchem CID for the specific compound.
        First tries a manually-compiled lookup table because chemists can't spell.
        Then attempts to find the CID by searching the extracted alias by name.
        Finally attempts to find the CID by searching the extracted alias as a formula.
        :param alias: a string representation identified by ChemicalTagger as a chemical.
        :return: None
        """
        assert isinstance(alias, str), f'I can only process chemical names as strings! {alias} is a {type(alias)}'

        try:
            pubchem_search = pcp.get_compounds(alias, 'name')
            self.pubchem_cid = pubchem_search[0].cid
            if len(pubchem_search) > 1:
                logging.warning(f'Multiple possible compounds for {alias} identified:\n{pubchem_search}')
        except IndexError:
            logging.warning(f'Following chemical is unidentified:\n {alias}')
            raise ChemicalNotFoundError
        except URLError as e:
            logging.warning(f'URL error when searching for compound {alias}\nFull text:\n{e}')
            raise ChemicalNotFoundError
        
        except pcp.PubChemHTTPError as e:
            sleep(5)
            try:
                pubchem_search = pcp.get_compounds(alias, 'name')
                self.pubchem_cid = pubchem_search[0].cid
            except:
                raise ChemicalNotFoundError
        except BadRequestError:
            raise ChemicalNotFoundError

        
    def find_pubchem_name_memoized(self, alias: str):
        if self.identifier_cache_location:
            try:
                with open(self.identifier_cache_location, 'r', encoding='utf-8') as f:
                    cache = json.loads(f.read())
            except (IOError, ValueError):
                logging.warn(f'No cache file found at {self.identifier_cache_location}')
                cache = {}
            
      
            try:
                self.pubchem_cid = cache[alias][0]
                logging.debug(f'identity cache at {self.identifier_cache_location} successfully used for {alias}')
            except KeyError:
                try: 
                    self.find_pubchem_name(alias)
                    cache[alias] = [self.pubchem_cid]
                    with open(self.identifier_cache_location, 'w', encoding='utf-8') as f:
                        json.dump(cache, f, indent=4)
                except ChemicalNotFoundError: 
                    cache[alias] = []
                    with open(self.identifier_cache_location, 'w', encoding='utf-8') as f:
                        json.dump(cache, f, indent=4)
                    raise ChemicalNotFoundError
            except IndexError:
                raise ChemicalNotFoundError


                
        else:
             logging.warn('No identity cache file specified!')
             self.find_pubchem_name(alias)
        
    def extract_key_properties_from_pubchem(self):
        """
        A helper function for pulling come basic information about a chemical from the pubchem database, given a CID.
        :return: None
        """
        pubchem_object = pcp.Compound.from_cid(self.pubchem_cid)
        self.name = pubchem_object.iupac_name
        if not isinstance(self.name, str):
            try:
                self.name = pubchem_object.synonyms[0]
            except IndexError:
                self.name = ''
            logging.debug(self.name)
        self.mol_wt = pubchem_object.molecular_weight
        self.mol_formula = pubchem_object.molecular_formula
        self.smiles = pubchem_object.isomeric_smiles
        self.inchi = pubchem_object.inchi

    def extract_key_properties_from_pubchem_memoized(self):
        if self.property_cache_location:
            try:
                with open(self.property_cache_location, 'r', encoding='utf-8') as f:
                    cache = json.loads(f.read())
            except (IOError, ValueError):
                logging.warn(f'No property cache file found at {self.property_cache_location}')
                cache = {}
            
            def update_cache(cache):
                cache[self.pubchem_cid] = {
                    'name': self.name,
                    'mol_wt': self.mol_wt,
                    'mol_formula': self.mol_formula,
                    'smiles': self.smiles,
                    'inchi': self.inchi
                }
      
            try:
                cache[self.pubchem_cid]
            except KeyError:
                self.extract_key_properties_from_pubchem()
                update_cache(cache)
                with open(self.property_cache_location, 'w', encoding='utf-8') as f:
                    json.dump(cache, f, indent=4)
            else:
                try:
                    self.name = cache[self.pubchem_cid]['name']
                    self.mol_wt = cache[self.pubchem_cid]['mol_wt']
                    self.mol_formula = cache[self.pubchem_cid]['mol_formula']
                    self.smiles = cache[self.pubchem_cid]['smiles']
                    self.inchi = cache[self.pubchem_cid]['inchi']
                    logging.debug(f'property cache at {self.property_cache_location} successfully used for {self.pubchem_cid}')
                except KeyError:
                    logging.warn(f'Invalid entry in cached properties for {self.pubchem_cid}! Defaulting to pubchem database search!')
                    self.extract_key_properties_from_pubchem()
                    update_cache(cache)
        else:
            logging.warn('No property cache specified!')
            self.extract_key_properties_from_pubchem()

    def find_chedl(self):
        """
        Finds the CAS identifier from a given INCHI code, for looking up the CHEDL database
        :return: None
        """
        try:
            working_metadata = chemicals.search_chemical(self.inchi)
            self.CAS = working_metadata.CASs
        except ValueError:
            pass

    def density_from_chedl(self, temp: float =298):
        """
        Estimates the density of a pure component for volume -> mole conversions.
        Uses the COSTALD algorithm for largely arbitrary reasons.
        :param temp: the temperature of the fluid, as a float (in K)
        :return: None
        """
        try:
            self.CAS
        except AttributeError:
            raise
        Tc = chemicals.critical.Tc(self.CAS)
        Pc = chemicals.critical.Pc(self.CAS)
        Zc = chemicals.critical.Zc(self.CAS)
        Vc = chemicals.critical.Vc(self.CAS)
        omega = chemicals.acentric.omega(self.CAS)

        COSTALD_density = chemicals.Vm_to_rho(chemicals.volume.COSTALD(T=temp, Tc=Tc, Vc=Vc, omega=omega),
                                              float(self.mol_wt))
        # print(f'COSTALD density = {COSTALD_density}')
        rackett_density = chemicals.Vm_to_rho(chemicals.volume.Rackett(T=temp, Tc=Tc, Pc=Pc, Zc=Zc),
                                              float(self.mol_wt))
        # print(f'rackett density = {rackett_density}')
        Yen_woods_density = chemicals.Vm_to_rho(chemicals.volume.Yen_Woods_saturation(T=temp, Tc=Tc, Vc=Vc, Zc=Zc),
                                                float(self.mol_wt))
        # print(f'yen_woods density = {Yen_woods_density}')

        self.density = Q_(float(COSTALD_density), ureg.g / ureg.liter).to_base_units()

    def to_json(self, cache_location='./'):
        try: 
            self.mol_wt
        except NameError:
            self.extract_key_properties_from_pubchem_memoized()
        output = {
        'name': self.pubchem_cid,
        'mol_wt': self.mol_wt,
        'mol_formula': self.mol_formula,
        'smiles': self.smiles,
        'inchi': self.inchi,
        }
        with open(f'{cache_location}pubchem_{self.pubchem_cid}.json', 'w') as f:
            f.write(json.dumps(output, indent=4))
    
    # TODO: write an explicit load from cache function (why?)
    # @staticmethod
    # def from_json(pubchem_cid, property_cache_location):
    #     with open(property_cache_location, 'r') as f:
    #         input = json.loads(f.read())
        

        
class ChemicalEntity(ChemicalSpecies):
    """
    A class containing primarily for the purposes fo relating a chemical identity to its amount.
    Contains functions to convert units into m0oles, for consistent analysis.

    Key functions:
    : calc_moles: unit converts an amount + unit into moles. Needs some refactoring
    : calc_moles_from_df_row: applies calc_moles to the row of a Pandas dataframe
    : calc_moles_from_string: applies calc_moles to a string
    : to_series: Outputs a well-formatted series for a BillOfMaterials

    TODO: refactor calc_moles workflow
    """



    def __init__(
            self,
            identifier: str =None,
            identifier_cache_location = None,
            property_cache_location = None
            ):
        super().__init__(
            alias=identifier, 
            identifier_cache_location=identifier_cache_location, 
            property_cache_location=property_cache_location
            )
        logging.info(f'Creating chemical entry for {self.pubchem_cid}: \n{str(self)}: ')

    def __str__(self):
        output = ''
        output += super().__str__()
        return output


    string_conversion_to_moles = {
        'mass': lambda quant, mol_wt, density: quant / mol_wt,
        'volume': lambda quant, mol_wt, density: quant * density / mol_wt,
        'other_amount': lambda quant, mol_wt, density: quant,
        # this will likely need fixing soon as other_amount is an ambigious term, but we'll see
        'conc': lambda quant, mol_wt, density: np.nan,
        # 'percent': lambda quant, mol_wt, density: np.nan
    }

    pint_conversion_to_moles = {
        UnitsContainer({'[mass]': 1}): string_conversion_to_moles['mass'],
        UnitsContainer({'[length]': 3}): string_conversion_to_moles['volume'],
        UnitsContainer({'[substance]': 1}): string_conversion_to_moles['other_amount']
    }



    def calc_moles(self, row_data):
        """
        Does a few too many things at once. Work flow should go like:
        1. Figure out the input units
        2. Calculate as moles
        3. Sum with any other values present

        :param row_data:
        :return:
        """
        logging.debug(f"working on the following row:\n{row_data}")

        # Catch incorrect data types in cell, return nan if they're not compatible
        try:
            if isinstance(row_data['Units used'], str):
                units_info = [row_data['Units used']]
            else:
                units_info = row_data['Units used']
            logging.debug(units_info)
        except (TypeError, KeyError):
            # logging.debug(row_data)
            self.moles = np.nan
            return

        # Create working total of moles and collect chemical molecular weight
        working_compound_total = Q_(0, ureg.mol)
        mol_wt = Q_(float(self.mol_wt), ureg.g / ureg.mol).to_base_units()

        # Confirm that units_info is iterable (i.e. contains a list of unit types to convert)
        try:
            iter(units_info)
        except TypeError:
            logging.debug(f'{units_info} is a {type(units_info)}, assuming this means no moles to calculate')
            self.moles = np.nan
            return

        # Create a Counter object of units to convert and dictionary of their strings
        units_used = Counter(units_info)

        # Create a dictionary of "unit type": [units as strings] to loop through unit-by-unit.
        to_convert = {}
        for k, v in units_used.items():
            if isinstance(row_data[k], str):
                working_values = [row_data[k]]
            else:
                working_values = row_data[k]
            logging.debug(f'{working_values}, {len(working_values)}')

            logging.debug(f'{k},{v}, {len(row_data[k])}, {row_data[k]}')
            if len(working_values) == v:
                to_convert[k] = working_values
            else:
                logging.warning(
                    f"There's a mismatch in the number of {k} units found and reported.\nThis entry will not be calculated just now.")
                continue

        logging.info(f"Calculating density now for {self.name}!")
        try:
            super().density_from_chedl()
            logging.debug(f'Density = {self.density}, {type(self.density)}')
            density = self.density
        except:
            density = Q_(np.nan, ureg.g / ureg.liter).to_base_units()
        logging.info(f'Calculated density for {self.name}: {density}')

        if self.pubchem_cid == 14923:  # special case for aqueous ammonia, because it's a nuisance
            density = Q_(900, ureg.g / ureg.liter).to_base_units()

        for k, v in to_convert.items():
            working_unit_total = 0
            try:
                for y in v:
                    logging.debug(y)
                    working_unit_total += ureg(y)
            except (UndefinedUnitError, TypeError, KeyError, DefinitionSyntaxError):
                logging.error("Error encountered converting raw info to pint unit:\n----------------------")
                logging.error(f'amount in question: {row_data}')
                break
            logging.debug(f'{working_unit_total}, {type(working_unit_total)}')
            try:
                working_compound_total += self.string_conversion_to_moles[k](working_unit_total,
                                                                             mol_wt, density).to(ureg.mol)
            except (TypeError, AttributeError, UndefinedUnitError):
                logging.error("Error encountered during calcualtions:\n----------------------")
                logging.error(f'amount in question: {y}')
                logging.error(f'Working unit total: {working_unit_total}')
                logging.error(f'Molecular weight:  {mol_wt}')
                logging.error(f'Density: {density}')

                raise
            except KeyError:
                working_unit_total = Q_(0, ureg.mol)

        logging.debug(row_data['name'], ': ', working_compound_total.to(ureg.mol))
        self.moles = round(working_compound_total.to(ureg.mol).magnitude, 6)
        logging.info(f'Moles calculated for {self.name}: {self.moles}')

    def calc_moles_from_df_row(self, row_data):
        self.calc_moles(row_data)

    def calc_moles_from_string(self, input_string):
        input_value = Q_(input_string)
        assert input_value.dimensionality in list(self.pint_conversion_to_moles.keys())

        mol_wt = Q_(float(self.mol_wt), ureg.g / ureg.mol).to_base_units()

        logging.info(f"Calculating density now for {self.name}!")
        try:
            super().density_from_chedl()
            logging.debug(f'Density = {self.density}, {type(self.density)}')
            density = self.density
        except:
            density = Q_(np.nan, ureg.g / ureg.liter).to_base_units()

        if self.pubchem_cid == 14923:  # special case for aqueous ammonia, because it's a nuisance
            density = Q_(900, ureg.g / ureg.liter).to_base_units()

        logging.debug(input_value, mol_wt, density)
        moles = self.pint_conversion_to_moles[input_value.dimensionality](input_value,
                                                                          mol_wt,
                                                                          density)

        self.moles = round(moles.to(ureg.mol).magnitude, 6)

    # TODO: create a calc moles from string interpreter Create a general function - if it's a string, attempt to
    #  infer units, if it's a dataframe row, attempt to infer from column info etc.

    def to_series(self):
        return pd.Series({'name': self.name, 'pubchem_id': int(self.pubchem_cid), 'moles': self.moles})

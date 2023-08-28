
"""
A module containing the tools to analyse synthesis Conditions and produce useful data structures for further analysis.

Author: Joe Manning (@jrhmanning, joseph.manning@manchester.ac.uk)
Date: Nov 2022

Classes:
Conditions - A class containing methods to build dataframes of reaction conditions

Exceptions:
None
"""
from __future__ import annotations
import pandas as pd
import re
import logging
import numpy as np
from pint import UndefinedUnitError, DimensionalityError, UnitRegistry, Unit

ureg = UnitRegistry()
Q_ = ureg.Quantity


class Conditions:
    """
    A class containing information about conditions in (part of) a synthesis sequence.
    Contains information about time, temperature, and space to expand e.g. into pH, pressure

    """

    def __init__(self, conditions: pd.DataFrame):
        self.raw_conditions = conditions
        self.extract_times_and_temps()
        self.extract_temps_only()
        self.extract_times_only()
    
    def extract_times_only(self):
        output = pd.DataFrame(columns=['step number','time'])
        try:
            input = self.raw_conditions[['step number','time']]
        except KeyError:
            input = self.raw_conditions['time']
        
        output_counter = 0
        for _,v in input.iterrows():
            if (isinstance(v['time'], str)) or (isinstance(v['time'], float)) and (v['time'] is not None):
                try:
                    working_time = eval(v['time'])
                except (SyntaxError, TypeError, NameError):
                    working_time = v['time']
            elif isinstance(v['time'], list):
                # print('Time is a list!')
                working_time = v['time']
            elif v['time'] is None:
                working_time = np.nan
            else:
                print(type(v['time']))
                raise SyntaxError

            try: 
                working_time
            except NameError: 
                raise(f"No temperature value found! Raw value is {v['time']}, ({type(v['time'])})")

            try:
                for x in working_time:
                    output.loc[output_counter, 'time'] = x
                    try:
                        output.loc[output_counter, 'step number'] = v['step number']
                    except: 
                        pass
                    output_counter += 1
            except TypeError:
                output.loc[output_counter, 'time'] = working_time
                output_counter += 1

        output['Time (min)'] = output['time'].apply(self.parse_time_to_minutes)
        self.times = output

    def extract_temps_only(self):
        output = pd.DataFrame(columns=['step number','temp'])
        try:
            input = self.raw_conditions[['step number','temp']]
        except KeyError:
            input = self.raw_conditions['temp']
        
        output_counter = 0
        for _,v in input.iterrows():
            if (isinstance(v['temp'], str)) or (isinstance(v['temp'], float)) and (v['temp'] is not None):
                try:
                    working_time = eval(v['temp'])
                except (SyntaxError, TypeError, NameError):
                    working_time = v['temp']
            elif isinstance(v['temp'], list):
                # print('Time is a list!')
                working_time = v['temp']
            elif v['temp'] is None:
                working_time = np.nan
            else:
                print(type(v['temp']))
                raise SyntaxError

            try: 
                working_time
            except NameError: 
                raise(f"No temperature value found! Raw value is {v['temp']}, ({type(v['temp'])})")

            try:
                for x in working_time:
                    output.loc[output_counter, 'temp'] = x
                    try:
                        output.loc[output_counter, 'step number'] = v['step number']
                    except: 
                        pass
                    output_counter += 1
            except TypeError:
                output.loc[output_counter, 'temp'] = working_time
                output_counter += 1

        output['T (K)'] = output['temp'].apply(self.parse_temperature_to_kelvin)
        self.temps = output            

    def extract_times_and_temps(self):
        """
        Produces a dataframe of quoted times and temperatures for a whole sequence in standardised units.
        TODO: add functionality to work for a partial sequence too.
        :return: None
        """
        try:
            output = self.raw_conditions[['step number','time', 'temp']]
        except KeyError:
            output = self.raw_conditions[['time', 'temp']]

        output.loc[:,'T (K)'] = np.nan
        output.loc[:,'T (K)'] = output.loc[:,'T (K)'].astype('O')
        output.loc[:,'Time (min)'] = np.nan
        output.loc[:,'Time (min)'] = output.loc[:,'Time (min)'].astype('O')


        for c, x in output.iterrows():
            logging.debug(c, x['temp'], x['time'])
                
            # Here I want to preprocess everything so that this function can handle strings as well as lists
            # So how would I?
            if (isinstance(x['temp'], str)) or (isinstance(x['temp'], float)) and (x['temp'] is not None):
                try:
                    working_temp = eval(x['temp'])
                except (SyntaxError, TypeError, NameError):
                    working_temp = [x['temp']]
            elif isinstance(x['temp'], list):
                working_temp = x['temp']
            elif x['temp'] is None:
                working_temp = []
            else:
                raise SyntaxError

            try: 
                working_temp
            except NameError: 
                raise(f"No temperature value found! Raw value is {x['temp']}, ({type(x['temp'])})")

            if (isinstance(x['time'], str)) or (isinstance(x['time'], float)) and (x['time'] is not None):
                try:
                    working_time = eval(x['time'])
                except (SyntaxError, TypeError, NameError):
                    working_time = [x['time']]
            elif isinstance(x['time'], list):
                # print('Time is a list!')
                working_time = x['time']
            elif x['time'] is None:
                working_time = []
            else:
                print(type(x['time']))
                raise SyntaxError

            try: working_time
            except NameError: raise(f"No temperature value found! Raw value is {x['time']}, ({type(x['time'])})")

            try:
                placeholder = [self.parse_temperature_to_kelvin(y) for y in working_temp]
                if len(placeholder)>0:
                    output.at[c, 'T (K)'] = placeholder
            except (IndexError, TypeError):
                pass

            try:
                placeholder = [self.parse_time_to_minutes(y) for y in working_time]
                if len(placeholder)>0:
                    output.at[c, 'Time (min)'] = placeholder
            except (IndexError, TypeError):
                pass


        self.time_temp = output

    temp_dict = {
        'C': ureg.degC,
        '°C': ureg.degC,
        'kelvin': ureg.kelvin,
        'Kelvin': ureg.kelvin,
        'k': ureg.kelvin,
        'K': ureg.kelvin
    }

    def parse_temperature_to_kelvin(self, df_temp_string: str) -> float:
        """
        Calculates the numerical temperature from a string of the temperature value + unit
        Performs some regex substitutions to recognise common phrases like "room temperature"
        :param df_temp_string: a string
        :return: The temperature, in kelvin
        """
        working = re.sub(r'\b(at|of|in|to|for)\b', '', df_temp_string).strip()
        working = re.sub(r'\b([0-9]+)-([0-9]+)\b', r'\1', working).strip()

        working = re.sub(r'[^\x00-\x7F]C', 'degC', working)
        working = re.sub('oC', 'degC', working)
        working = re.sub(r'\bC', 'degC', working)

        try:
            output = Q_(float(working.split()[-2]), self.temp_dict[working.split()[-1]])
        except (ValueError, UndefinedUnitError, IndexError) as e:
            if bool(re.search(r'(Room|room|ambient|indoor)\Wtemperature|RT|(A|a)mbient', working)):
                output = Q_(25, ureg.parse_units('degC'))
            else:
                logging.warn(f'Failed for: "{working}"')
                logging.warn(f'\t {type(e).__name__} error raised: {e}')
                return np.nan
        except KeyError as e:
            try:
                output =  Q_(float(working.split()[-2]), working.split()[-1])
            except Exception as e:
                print(f'Temperature extraction failed for: "{working}", original was "{df_temp_string}"')
                print(type(e).__name__, e)
                return np.nan
        except Exception as e:
            print(f'Temperature extraction failed for: "{working}", original was "{df_temp_string}"')
            print(type(e).__name__, e)
            return np.nan
        # print(output)
        try:
            return output.to_base_units().magnitude
        except DimensionalityError:
            return np.nan

    def parse_time_to_minutes(self, df_time_string: str) -> float:
        """
        Identifies a numerical time from a string of time quantity + unit.
        Performs regex to replace most common values for pint-parseable alternatives.
        :param df_time_string: a timephrase string
        :return: the amount of minutes used as a float
        """
        working = re.sub(r'\b(at|of|in|to|for|(A|a)fter)\b', '', df_time_string).strip()
        working = re.sub(r'~|∼', '', working).strip()
        working = re.sub(r'\b(h)\b', 'hour', working).strip()
        working = re.sub(r'\b(d|D)\b', 'day', working).strip()
        working = re.sub(r'\b(P|p)eriod\b', '', working).strip()


        working = re.sub(r'\b(an|a|one)\W(hour|day)\b', r'1 \2', working).strip()
        working = re.sub(r'\b(an|a|one)\W(night)\b', r'18 hours', working).strip()
        working = re.sub(r'\b(one)\b', r'1 ', working).strip()
        working = re.sub(r'\b(two)\b', r'2 ', working).strip()
        working = re.sub(r'\b(three)\b', r'3 ', working).strip()
        working = re.sub(r'\b(four)\b', r'4 ', working).strip()
        working = re.sub(r'\b(five)\b', r'5 ', working).strip()
        working = re.sub(r'\b(six)\b', r'6 ', working).strip()
        working = re.sub(r'\b(seven)\b', r'7 ', working).strip()
        working = re.sub(r'\b(eight)\b', r'8 ', working).strip()
        working = re.sub(r'\b(nine)\b', r'9 ', working).strip()
        working = re.sub(r'\b(ten)\b', r'10 ', working).strip()
        working = re.sub(r'\b(eleven)\b', r'11 ', working).strip()
        working = re.sub(r'\b(twelve)\b', r'12 ', working).strip()

        working = re.sub(r'\b([0-9.]+)(-|~)([0-9.]+)\b', r'\1', working).strip()

        working = re.sub(r'\b(mins)\b', r'minutes ', working).strip()

        if bool(re.search('(O|o)ver.?night|(O|o)ne night', working)):
            output = Q_(18, ureg.hour)

        else:
            try:
                output = Q_(float(working.split()[-2]), working.split()[-1])
            except (ValueError, KeyError, UndefinedUnitError, IndexError, TypeError) as e:
                print(f'Time extraction failed for: "{working}", original was "{df_time_string}"')
                print(type(e).__name__ ,e)
                return np.nan
            except DimensionalityError:
                return np.nan
        try:
            return output.to('minutes').magnitude
        except DimensionalityError:
            return np.nan

    def calculate_similarity(self, other: Conditions) -> float:
        return (self._temp_similarity(other)+self._time_similarity(other))/2
        
    def _temp_similarity(self, other: Conditions) -> float:
        """
        Simple similarity calculatoin based on (relative) number of dissimilar temperature values quoted duirng a synthesis.
        Liable to give weird values if performed on an entire sequence at once e.g. if room temperature is quoted in two separate places.
        However should (hopefully) be reliable when applied to a small enough subsection of a sequence.
        """
        clean_self_set = set(self.temps['T (K)'][self.temps['T (K)'].notna()].to_list())
        clean_other_set = set(other.temps['T (K)'][other.temps['T (K)'].notna()].to_list())

        similar_temperatures = list(clean_self_set.intersection(clean_other_set))
        all_temperatures = list(clean_self_set.union(clean_other_set))

        dissimilar_temperatures_self = clean_self_set.difference(clean_other_set)
        dissimilar_temperatures_other = clean_other_set.difference(clean_self_set)


        try:
            return 100*(len(similar_temperatures)/len(all_temperatures))
        except ZeroDivisionError:
            return 0

    def _time_similarity(self, other: Conditions) -> float:
        """
        Simple similairty calculation based on relative difference in total times quoted. 
        To provide more detail/granularity this method requires conisderaiton of multiple subsections of a sequence.
        """
        current_total_time = self.times['Time (min)'][self.times['Time (min)'].notna()].sum()
        other_total_time = other.times['Time (min)'][other.times['Time (min)'].notna()].sum()

        try:
            return 100*(1-(abs(current_total_time-other_total_time)/((current_total_time+other_total_time)/2)))
        except ZeroDivisionError:
            return 0
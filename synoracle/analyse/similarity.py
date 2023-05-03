"""
A few methods to determine the completeness of a unit operation based on its contents.
Aims to provide a score between 0 and 1 fo a 

Author: Joe Manning (@jrhmanning, joseph.manning@manchester.ac.uk)
Date: Jan 2022

Key methods:

"""

import pandas as pd
import logging
from ..interpret.ingredients import ChemicalList, BillOfMaterials, EmptyChemicalList
from ..interpret.conditions import Conditions


def calculate_completeness(df_row) -> float:
    """
    A function to calculate the completeness of a synthesis sequence "block" based on flexible criteria.
    Takes in a single row of a sequence dataframe containing "Step supertype", "new_chemicals", "time", "temperature" columns.
    Matches the step supertype against a specific set of criteria in a dictionary.
    These criteria are functions returning a number between 0 and 1 where 0 is wholly incomplete and 1 is wholly complete.
    For example an "add" step requires all chemicals to have a quantity to be meaningful. Therefore 0 = no chemical quantities and 1 = all chemical quantities.
    """

    suptertype_categories = {
        'add': calculate_bom_completeness,
        'react': calculate_condition_completeness,
        'remove': lambda x: (calculate_condition_completeness(x)+calculate_bom_completeness(x))/2
    }

    try:
        return suptertype_categories[df_row['Step supertype']](df_row)
    except KeyError:
        return 0


def calculate_bom_completeness(df_row) -> float: 
    """
    Calculates the completeness of a bill of materials inside a sequence.
    Generates the bill of mateials for th dataframe row, then calculates the fraction of identified chemical entities with 
    """
    try:
        chemlist = ChemicalList(pd.DataFrame(df_row['new_chemicals']))
    except EmptyChemicalList:
        return 0

    chemlist.prioritise_bom_df()
    bill = chemlist.produce_bill_of_mats().bill_of_materials
    # bill = BillOfMaterials(chemlist.chemical_list_reduced)
    try:
        completeness = len(bill['moles']>0)/len(bill)
    except ZeroDivisionError:
        completeness = 0

    return completeness 

def calculate_condition_completeness(df_row) -> float:
    """
    Calculates the completess of a block's time and temperature conditions

    """
    conditions = Conditions(pd.DataFrame(df_row).T)
    try:
        time_completeness = len(conditions.time_temp['Time (min)']>0)/len(conditions.time_temp)
    except (KeyError, ZeroDivisionError):
        time_completeness = 0
    try:
        temp_completeness = len(conditions.time_temp['T (K)']>0)/len(conditions.time_temp)
    except (KeyError, ZeroDivisionError):
        temp_completeness = 0

    completeness = (temp_completeness*time_completeness)/2

    return completeness

"""
A module containing classes for handling synthesis paragraphs as both text and XML documents (and converting between).
Makes use of ChemicalTagger to turn text paragraphs into hierarchical XML trees, then extracts them to pandas dataframes.

Author: Joe Manning (@jrhmanning, joseph.manning@manchester.ac.uk)
Date: Nov 2022

Classes:
SynParagraph

Exceptions:
InputFileContentError
InvalidInputError
"""
from lxml import etree
from lxml.etree import XMLSyntaxError, _Element
import logging
import numpy as np
from pathlib import Path
import pandas as pd
from copy import deepcopy
from typing import Union
from itertools import chain
import re
from typing import List, Tuple

# create logger
logger = logging.getLogger('simple_example.txt')
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)

class InputFileContentError(Exception): pass

class InvalidInputError(Exception): pass


class SynParagraph:
    '''
    An object for taking in a raw text paragraph and producing a structured synthesis sequence.

    TODO: add methods to handle identified ratios and also punctuation-delimited lists of chemicals

    Key methods:
    : regex_preprocess: performs some quality-of-life improvements to help ChemicalTagger correctly split tokens
    : apply_chem_tagger: runs the ChemicalTagger Java executable externally to transform the data into hierarchical xml
    : load_xml: Loads in the resulting XML as an element tree for sequential analysis
    : xml_para_annotate: Some fun text annotation to demonstrate the ChemicalTagger actions
    : find_chemical_name: creates a list of names used for a single chemical compound within an XML tag
    : find_chemical_quantity: creates a string of quantity information for a specific unit type within an XMl tag
    : find_chemicals: Associates chemical names and quantities within an XML tag
    : parse_actionphrase: Identifies all chemicals and conditions from within a specific Action XML tag
    : process_actionphrases: Iterates through all action tags in a sentence to create a sequence of verbs and features
    : extract_sequence: Performs the above for a full paragraph, creating a Dataframe of a Sequence

    '''

    def __init__(self, paper_identifier: str, source_directory: Union[str, Path] = Path('./')):
        """
        Instantiates the object and concerts a text document to XML (if needed)

        :param paper_identifier: a unique string pointing to the synthesis paragraph as a text file
        :param source_directory: a string or path pointing to the directory where your input file is
        """
        self.source_directory = Path(source_directory)
        print(source_directory)
        self.paper_indentifier = paper_identifier
        self.source_paragraph = self.source_directory / (paper_identifier + '.txt')
        #self.regex_preprocess()
        self.load_xml()
        self.extract_sequence()

    def regex_preprocess(self):
        """
        Replaces all instances of specific tokens within a text synthesis with cleaned versions
        TODO: make the outputfile slightly more sensible like x.##.cleaned.txt
        :return: None
        """
        with open(self.source_paragraph, 'r', encoding='utf-8') as f:
            rawtext = f.read()
            rawtext2 = re.sub(r'([0-9]+)(K)\b', r'\1 \2', rawtext)
        with open('placeholder.txt', 'w', encoding='utf-8') as f2:
            f2.write(rawtext2)

    def apply_chem_tagger(self, chemtagger_dir: Union[str, Path] = './',
                          chemtagger_exec: str = 'chemicalTagger-1.6-SNAPSHOT-jar-with-dependencies-file.jar') -> Path:
        """
        Applied ChemicalTagger to a specific paragraph, if an xml with the same name doesn't yet exist
        :param chemtagger_dir: The chemtagger executable location
        :param chemtagger_exec: Name of the chemtagger executable, in case you changed yours
        :return: the directory path for the xml file generated
        """
        logging.debug(f"Applying chem tagger on {self.paper_indentifier}")

        paragraph = self.source_paragraph
        function_output = self.source_directory / (self.paper_indentifier + '.xml')

        try:
            assert function_output.is_file(), paragraph
        except AssertionError:
            raise NotImplementedError
        return function_output

    def load_xml(self):
        """
        Loads an XML file into memory as an ElementTree
        :return: None
        """
        xml_filename = self.apply_chem_tagger()
        if not xml_filename.is_file():
            raise InvalidInputError(f"Cannot find extracted xml actions for paper {self.paper_indentifier}")
        with open(xml_filename, 'rb') as f:
            raw = f.read()
        try:
            self.working_xml = etree.fromstring(raw)
        except XMLSyntaxError as e:
            logging.error('Cannot read extracted xml actions for paper {0}'.format(self.paper_indentifier))
            raise InputFileContentError
        logging.info('XML loaded in from {0}'.format(xml_filename.parts[-1]))

    def _text_annotate(self, text: str, start_char_list: list, end_char_list: list, texttype: list = ['bold']) -> str:
        """
        Provides unicode string formatting for printing strings with the ChemicalTagger data annotations highlighted
        :param text: The raw text string to modify
        :param start_char_list: The position of characters to start the given formatting at
        :param end_char_list: The position of characters to end the formatting at
        :param texttype: The formatting type required
        :return: a unicode-modified version of the input string
        """
        format_chars = {
            'purple': '95',
            'cyan': '96',
            'darkcyan': '36',
            'blue': '94',
            'green': '92',
            'yellow': '93',
            'red': '91',
            'bold': '1',
            'underline': '4',
            'italic': '3',
            'strikethrough': '9',
            'reverse': '7',
            'end': '0'
        }
        if isinstance(texttype, str):
            texttype = [texttype]
        texttype = [x.lower() for x in texttype]
        assert all([format_chars[x] for x in texttype]), f'Invalid text formatting choice {texttype}.'
        start_points = sorted(start_char_list, reverse=True)
        end_points = sorted(end_char_list, reverse=True)
        for start_char, end_char in zip(start_points, end_points):
            text = text[:start_char] + '\033[' + ';'.join([format_chars[x] for x in texttype]) + 'm' + text[
                                                                                                       start_char:end_char] + '\033[' + \
                   format_chars['end'] + 'm' + text[end_char:]
        return text

    def xml_para_annotate(self, raw_xml: etree.ElementTree) -> str:
        """
        Produces an annotated string exemplifying the XML tags included in the XMl sequence
        :param raw_xml: The XML element tree to be annotated
        :return: a string of annotated text
        """
        output = ''
        tree = etree.ElementTree(raw_xml)
        colour_typing = {
            'Synthesize': 'purple',
            'Dissolve': 'cyan',
            'Add': 'red',
            'Dry': 'darkcyan',
            'Stir': 'blue',
            'Heat': 'green',
            'Wash': 'yellow',
        }
        for x in tree.iter():
            if x.text is not None:
                parents = tree.getpath(x)
                annotations = []
                if 'OSCAR-CM' in parents:
                    annotations.append('underline')
                if 'QUANTITY' in parents:
                    annotations.append('bold')
                if 'ActionPhrase' in parents:

                    above = x.getparent()
                    while above.tag not in ['ActionPhrase', None]:
                        above = above.getparent()
                    try:
                        annotations.append(colour_typing[above.attrib['type']])
                    except KeyError:
                        annotations.append('purple')
                working = self._text_annotate(x.text + ' ', [0], [-1], annotations)

                output += working

        return output

    def find_chemical_name(self, xml: _Element) -> Tuple[list, str]:
        """
        Produces the chemical name from an XML COMPOUND tag. based on OSCAR named entity recognition.
        Currently just joins all the tokens together for each OSCARCM tag.
        :param xml: an XML element describing containing chemical name(s)
        :return: a tuple of a list of names and the longest one found
        """
        assert type(xml) == _Element
        from itertools import chain

        # process the result of xml.iter(tag="MOLECULE") or xml.iter(tag="UNNAMEDMOLECULE")
        aliases = []
        for names in chain(xml.iter(tag="OSCARCM"), xml.iter(tag="REFERENCETOCOMPOUND")):
            aliases.append(' '.join(list(names.itertext())))
        try:
            mol_name = max(aliases) # I don't kno why I used max here. Maximum name length perhaps?
        except ValueError:
            mol_name = 'unknown'

        return aliases, mol_name

    def find_chemical_quantity(self, xml: _Element, tag: str) -> str:
        """
        Outputs the chemical quantity mentioned in an XML tag as a string for later analysis.
        Currently can only accept one quantity tg of each type mentioned.
        Should this functionality be extended?
        :param xml: an XML element containing chemical quantity information
        :param tag: Type of physical quantity to be analsyed
        :return: a string of the chemical quantity mentioned
        """
        assert type(xml) == _Element

        # process the result of xml.iter(tag="QUANTITY")
        # tag should be one of 'MASS', 'AMOUNT', 'VOLUME', 'PERCENT', 'MOLAR'
        if len(list(xml.iter(tag=tag))) > 1:
            logging.warning(f"Warning! Multiple ({len(list(xml.iter(tag=tag)))}) {tag} tags found!")
            return np.nan
        if len(list(xml.iter(tag=tag))) == 0:
            return np.nan
        return ' '.join([' '.join(list(x.itertext())) for x in xml.iter(tag=tag)])

    def find_chemicals(self, xml: _Element) -> List[dict]:
        """
        Iterates though an ActionPhrase to find all the mentions of chemicals and their quantities therein.
        :param xml: An ActionPhrase tag potentially containing chemicals
        :return: a list of chemical dictionaries containing information on name, aliases, and amounts of various types.
        """
        assert type(xml) == _Element

        # process the result of xml.iter(tag="ActionPhrase") for chemicals
        from itertools import chain
        outputs = []
        for chemical in chain(xml.iter(tag="MOLECULE"), xml.iter(tag="UNNAMEDMOLECULE")):
            aliases, mol_name = self.find_chemical_name(chemical)
            molec = {
                'name': mol_name,
                'mass': self.find_chemical_quantity(chemical, tag='MASS'),
                'other_amount': self.find_chemical_quantity(chemical, tag='AMOUNT'),
                'volume': self.find_chemical_quantity(chemical, tag='VOLUME'),
                'percent': self.find_chemical_quantity(chemical, tag='PERCENT'),
                'concentration': self.find_chemical_quantity(chemical, tag='MOLAR'),
                'aliases': aliases
            }
            outputs.append(molec)
        return outputs

    def parse_actionphrase(self, xml: _Element, counter: int) -> Tuple[dict, int]:
        """
        Processes a single action phrase for specific features and information of chemicals mentioned.
        :param xml: The ActionPhrase in question for processing
        :param counter: the sequential position of the actionphrase within the text
        :return: a dictionary of information about the actionphrase and an incremented counter index
        """
        assert type(xml) == _Element
        try:
            output_dict = {counter: {'name': xml.attrib['type'],
                                     'text': ' '.join(list(xml.itertext())),
                                     'new_chemicals': self.find_chemicals(xml),
                                     'temp': [' '.join(list(x.itertext())) for x in xml.iter(tag="TempPhrase")],
                                     'time': [' '.join(list(x.itertext())) for x in xml.iter(tag="TimePhrase")],
                                    'prepphrase': [' '.join(list(x.itertext())) for x in xml.iter(tag="PrepPhrase")],
                                    'apparatus': [' '.join(list(x.itertext())) for x in chain(xml.iter(tag="APPARATUS"), xml.iter(tag="VB-APPARATUS"))],

                                     'step number': counter
                                     }}
        except KeyError:
            output_dict = {counter: {'name': None,
                                     'text': ' '.join(list(xml.itertext())),
                                     'new_chemicals': self.find_chemicals(xml),
                                     'temp': [' '.join(list(x.itertext())) for x in xml.iter(tag="TempPhrase")],
                                     'time': [' '.join(list(x.itertext())) for x in xml.iter(tag="TimePhrase")],
                                    'prepphrase': [' '.join(list(x.itertext())) for x in xml.iter(tag="PrepPhrase")],
                                      'apparatus': [' '.join(list(x.itertext())) for x in chain(xml.iter(tag="APPARATUS"), xml.iter(tag="VB-APPARATUS"))],

                                     'step number': counter
                                     }}
        counter += 1
        return output_dict, counter

    def process_actionphrases(self, xml: _Element, counter: int=0) -> Tuple[dict, int]:
        """
        Performs parse_actionphrases on a full sentence, providing a sequential list of decribes actions.
        Operates on a depth-first approach for nested actions. I.E. for the case of "I poured a solution of X in Y",
        the output becomes 0: dissolve X in Y, 1: pour.
        :param xml: An XML element of a full sentence containing multiple possible ActionPhrases
        :param counter: the sequential position of the actionphrase within the text
        :return: a dictionary of actionphrases with appropraitely incremented action counter
        """
        logging.debug("Processing actionphrases in xml tag {0}, starting from action {1}".format(xml.tag, counter))
        assert type(xml) == _Element
        output = {}
        for action in xml.iter(tag='ActionPhrase'):
            logging.debug('Processing action of type {0}'.format(action.attrib))
            for child in action:
                for subaction in child.iter(tag='ActionPhrase'):
                    logging.debug('Processing subaction of type {0}'.format(subaction.attrib))
                    placeholder, counter = self.parse_actionphrase(subaction, counter)
                    output = {**output, **placeholder}
                    logging.debug('Subaction appended, current sequence is:')
                    logging.debug('\t {0}'.format([f'{k}: {v["name"]}' for k, v in output.items()]))
                    subaction.clear()
            placeholder, counter = self.parse_actionphrase(action, counter)
            output = {**output, **placeholder}
            logging.debug('Action appended, current sequence is:')
            logging.debug('\t {0}'.format([f'{k}: {v["name"]}' for k, v in output.items()]))
            action.clear()
        unassigned = {
            'chems': self.find_chemicals(xml),
            'times': [' '.join(list(x.itertext())) for x in chain(xml.iter(tag="TimePhrase"), xml.iter(tag="NN-TIME"))],
            'temps': [' '.join(list(x.itertext())) for x in xml.iter(tag="TempPhrase")],
            'prepphrase': [' '.join(list(x.itertext())) for x in xml.iter(tag="PrepPhrase")],
            'apparatus': [' '.join(list(x.itertext())) for x in chain(xml.iter(tag="APPARATUS"), xml.iter(tag="VB-APPARATUS"))]
        }
        if any([len(x) > 0 for x in unassigned.values()]):
            output[counter] = {
                'name': None,
                'text': ' '.join(list(xml.itertext())),
                'new_chemicals': unassigned['chems'],
                'temp': unassigned['temps'],
                'time': unassigned['times'],
                'prepphrase': unassigned['prepphrase'],
                'apparatus': unassigned['apparatus'],
                'step number': counter
            }
            counter += 1
        return output, counter

    def extract_sequence(self):
        """
        Performs process_actionphrases across all sentences within a paragraph, outputting as a pandas DataFrame
        :return: None
        """
        raw_sequence = {}
        counter = 0
        placeholder = deepcopy(self.working_xml)
        for x in placeholder.findall('Sentence'):
            output, counter = self.process_actionphrases(x, counter)
            raw_sequence = {**raw_sequence, **output}
        self.raw_synthesis = pd.DataFrame(raw_sequence).T

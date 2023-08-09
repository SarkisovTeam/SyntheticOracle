<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->




<!-- PROJECT LOGO -->

[//]: # (<br />)

[//]: # (<div align="center">)

[//]: # (  <a href="https://github.com/github_username/repo_name">)

[//]: # (    <img src="images/logo.png" alt="Logo" width="80" height="80">)

[//]: # (  </a>)

<h3 align="center">Synthetic Oracle</h3>
<div>
  <p align="center">
    A python-based tool for parsing and interpreting experimental synthesis reports.
    <br />
    <!-- <a href="https://github.com/SarkisovTeam/SyntheticOracle/doc"><strong>Explore the docs »</strong></a> -->
    <br />
    <a href="https://github.com/SarkisovTeam/SyntheticOracle/tree/main/worked%20example">View Worked Example</a> 
    <!-- <br /> -->
    <!-- <a href="https://github.com/SarkisovTeam/SyntheticOracle/demo">View Demo</a> -->
    ·
    <a href="https://github.com/SarkisovTeam/SyntheticOracle/issues">Report Bug</a>
    ·
    <a href="https://github.com/SarkisovTeam/SyntheticOracle/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

[//]: # ([![Product Name Screen Shot][product-screenshot]]&#40;https://example.com&#41;)

A Python-based tool for large-scale analysis and comparison of synthesis protocols. This module contains tools to:

1. Organise a sequence of synthesis actions from a well structured .json object (typically generated as an output from our [preprocessing workflow](https://github.com/SarkisovTeam/SynOracle-preprocessing)
2. Cross-reference chemical names against the PubChem database of chemical entities
3. Standardise all physical quantities into their SI units
5. Perform quantitative meta-analysis on the resultant corpus of synthesis protocols

### Built with
* [Pubchempy](https://pubchempy.readthedocs.io/en/latest/) [(Pypi)](https://pypi.org/project/PubChemPy/1.0/) [(Git)](https://github.com/mcs07/PubChemPy)
* [Pandas](https://pandas.pydata.org/) [(Pypi)](https://pypi.org/project/pandas/) [(Git)](https://github.com/pandas-dev/pandas)
* [Pint](https://pint.readthedocs.io/en/stable) [(Pypi)](https://pypi.org/project/Pint/) [(Git)](https://github.com/hgrecco/pint)


<!-- Here's a blank template to get started: To avoid retyping too much info. Do a search and replace with your text editor for the following: `github_username`, `repo_name`, `twitter_handle`, `linkedin_username`, `email_client`, `email`, `project_title`, `project_description` -->

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

This branch of the software is python native, requiring python 3.7+ environment (although there are issue swith python 3.10 as of 31-1-23). 
All required packages can be installed through pip with the following command:

   ```sh
   pip install -r requirements.txt
   ```

### Installation

1. Clone the repository
   ```sh
   git clone https://github.com/SarkisovTeam/SyntheticOracle.git
   ```
2. Install pip requirements (preferable in a fresh python environment)
   ```sh
   pip install -r requirements.txt
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

The software reads in raw synthesis sequences from JSON objects with the following tags titles:

```json
{
  "Step name": {
      "0": "Add"
   },
  "text": {
    "0": "Some raw text here"
  },
  "new_chemicals": {
    "0": [
      {
      "name": "Methanol", 
      "mass": "0.5 g", 
      "other_amount": "0.2 mmol", 
      "volume": "0.6 mL"
      }
    ]
  },
  "temp": {
    "0": ["at 300 oC"]
  },
  "time": {
    "0": ["for 24 hours"]
  }
}
```

Descriptions of how to produce these data structures are provided in the companion repository [here](https://github.com/SarkisovTeam/SynOracle-preprocessing). Once generated, the JSON data can be used to instantiate a `Sequence` object. Through the `extract_chemicals` class method, a `ChemicalList` object is created containing a set of records for each chemical mentioned in the sequence, which can be converted to a summarised `BillOfMaterials` using the `produce_bill_of_mats` method. Finally, the `extract_conditions` method of the `Sequence` object can produce a table of synthesis times and temperatures (a `Conditions` object). A worked example of the data workflow is shown in the <a href="https://github.com/SarkisovTeam/SyntheticOracle/example">example</a> folder.

<!-- 

To build a database of published syntheses for a given material, a six-step workflow is required, which will be briefly described here. Jupyter notebooks are provided in  with example code for ZIF-8, for your convenience.

1. Locate papers using `elsapy` and your self-defined keywords
2. Download a corpus of papers from each publisher identified
3. Extract plain text synthesis paragraphs from each paper downloaded
4. Process the paragraphs into hierarchical XML sequences, and extract sequential information to dataframes
5. Cross-reference extracted data against chemical databases to calculate quantities in standardised units
6. Summarise and analyse the identified sequences for further analysis

Each of these steps are discussed in detail in the documentation and demo. For convenience, minimum working examples for each step are provided below. 

### Finding papers


### Downloading papers
Generally speaking, each publisher has their own internal rules for text and data mining (TDM), including unique file formats for each individual paper as well as unique methods of accessing them. This document won’t provide an exhaustive guide (which can be found elsewhere) but will attempt to give an overview of the field using case studies from some of the largest chemical publishers. In general, there are three different strategies to perform TDM through a publisher – by interfacing with a human, with the publisher’s website directly, or with an online REST API, or. Each of these methods is contingent on your institution having a subscription with the publisher in question, which you check by trying to manually access the paper using your institutional login. 


### Identifying synthesis paragraphs

### Extracting hierarchical data from paragraphs

### Cross-referencing chemical information extracted

### Analysing data trends

-->

<!-- Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_
-->
<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ROADMAP -->
## Roadmap

### Short term (<1 month )

- [ ] Well-documented demo files (and/or jupyter notebooks) demonstrating code features and functionality
- [ ] Jupyter notebook examples to reproduce informaiton and figures from the associated manuscript.
- [ ] Pruning old and unused data

### Medium-term (< 12 months)
- [ ] Identification of MOFID to automatically identify chemical constraints on synthesis
- [ ] Quality metrics for synthesis parsing
- [ ] Generating material-by-material synthesis reports with summary data and details, à la David Fairen-Jiminez' [MOFexplorer](http://aam.ceb.cam.ac.uk/mofexplorer.html)
    - [ ] Define steps to build a plotly dashboard
    - [ ] Put them in here

See the [open issues](https://github.com/github_username/repo_name/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Joe Manning - [@jrhmanning](https://twitter.com/jrhmanning) - joseph.manning@manchester.ac.uk

<!-- Project Link: [https://github.com/github_username/repo_name](https://github.com/github_username/repo_name) -->

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
<!--
## Acknowledgments

* []()
* []()
* []()

-->
<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo_name.svg?style=for-the-badge
[contributors-url]: https://github.com/github_username/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo_name.svg?style=for-the-badge
[forks-url]: https://github.com/github_username/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo_name.svg?style=for-the-badge
[stars-url]: https://github.com/github_username/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo_name.svg?style=for-the-badge
[issues-url]: https://github.com/github_username/repo_name/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo_name.svg?style=for-the-badge
[license-url]: https://github.com/github_username/repo_name/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 

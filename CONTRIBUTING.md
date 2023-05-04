# Contributing to ElmerFEM

:+1::tada: Thank you for your interest to contribute to ElmerFEM! 🎉::+1:

The following is a set of guidelines for contributing to Elmer, on [GitHub](https://github.com/ElmerCSC/elmerfem). Without wanting to impose too much constraints on your contribution, please see this documents as a general guideline.

#### Table Of Contents

[Code of Conduct](#code-of-conduct)

[I just have a question to ask about Elmer!](#i-just-have-a-question)

[What should I know before I get started?](#what-should-i-know-before-i-get-started)
  * [Elmer code development philosophy](#elmer-code-development-philosophy)


[How Can I Contribute?](#how-can-i-contribute)
  * [Reporting Bugs](#reporting-bugs)
  * [Suggesting Enhancements](#suggesting-enhancements)
  * [Your First Code Contribution](#your-first-code-contribution)
  * [Pull Requests](#pull-requests)

[Styleguides](#styleguides)
  * [F90 style guide](#f90-style-guide)
  * [Git Commit Messages](#git-commit-messages)
  * [Documentation Styleguide](#documentation-styleguide)

[Additional Notes](#additional-notes)
  * [Issue and Pull Request Labels](#issue-and-pull-request-labels)

## Code of Conduct

This project and everyone participating in it is governed by the [ElmerFEM Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behaviour to the [ElmerFEM Administrator](mailto:elmeradmin@csc.fi).

## I just have a question!

> **Note:** Please refrain from opening an issue here at GitHub to ask a question.] 

We maintain an official [Elmer forum](https://www.elmerfem.org/forum/) - that's exactly the place were to put questions (or most likely already find the answers posted)

## What should I know before I get started?
### Elmer code development philosophy
We expect you to be familiar with the code development philosophy of Elmer. The up-to-date version of Elmer is in repository branch `devel`, which by name is a constantly changing code-base. Larger separate steps in development take place in separate feature branches that in the end shall be merged into `devel`. If  a new feature to an existing model or even a complete new model we expect a small test case that can be integrated into our `ctest` test suite that provides a benchmark for the new functionality. Any further new feature to the library and the modules of will then be evaluated against that suite of tests.

## How Can I Contribute?

### Reporting Bugs

This section guides you through submitting a bug report for ElmerFEM. Following these guidelines helps maintainers and the community understand your report :pencil:, reproduce the behaviour :computer: :computer:, and find related reports :mag_right:.

Before creating bug reports, please check [this list](#before-submitting-a-bug-report) as you might find out that you don't need to create one. When you are creating a bug report, please [include as many details as possible](#how-do-i-submit-a-good-bug-report).

> **Note:** If you find a **Closed** issue that seems like it is the same thing that you're experiencing, open a new issue and include a link to the original issue in the body of your new one.

#### Before Submitting A Bug Report

In order to avoid multiple reporting of bugs, please follow these steps before reporting one
* make sure you are on the correct branch (we only take reports from ´devel`) and that it is up-to-date
* check the current issues filed for this branch to whether the problem already has been mentioned
* make sure the fault of the code is not caused by errors in the inputs

#### How Do I Submit A (Good) Bug Report?

Bugs are tracked as [GitHub issues](https://guides.github.com/features/issues/). 

Explain the problem and include additional details to help maintainers reproduce the problem:

* **Use a clear and descriptive title** for the issue to identify the problem.
* **Describe the exact steps which reproduce the problem** in as many details as possible. When listing steps, **don't just say what you did, but explain how you did it**.
* **Provide specific examples to demonstrate the steps**. Include links to files or GitHub projects, or copy/pasteable snippets, which you use in those examples. If you're providing snippets in the issue, use [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).
* Best would be to **provide the (hopefully small) case that reproduces the behaviour** on your side.
* **Describe the behaviour you observed after following the steps** and point out what exactly is the problem with that behaviour.



Include details about your configuration and environment:

* **Which version of Elmer are you using?** You should see the version plus the date when compiled in the first output lines of your Elmer simulation,.
* **What's the name and version of the OS you are running Elmer on**?
* What compiler suite and what flags you used to compile Elmer?
* What hardware Elmer is run on?
* What external libraries are linked into your Elmer version, in particular the MPI version?


### Suggesting Enhancements

This section guides you through submitting an enhancement suggestion for Elmer, including completely new features and minor improvements to existing functionality. Following these guidelines helps maintainers and the community understand your suggestion :pencil: and find related suggestions :mag_right:.

Before creating enhancement suggestions, please check [this list](#before-submitting-an-enhancement-suggestion) as you might find out that you don't need to create one.

#### Before Submitting An Enhancement Suggestion

* **Check whether this enhancement has already been suggested**.
* Check, whether there is a dedicated branch of Elmer that already contains that enhancement or has it in development.

#### How Do I Submit A (Good) Enhancement Suggestion?

Enhancement suggestions are tracked as [GitHub issues](https://guides.github.com/features/issues/). 

* **Use a clear and descriptive title** for the issue to identify the suggestion.
* **Provide a step-by-step description of the suggested enhancement** in as many details as possible.
* **Specify the name and version of the OS you are using.**

### Your First Code Contribution



### Pull Requests

The process described here has several goals:

- Maintain ElmerFEM's quality
- Fix problems that are important to users
- Engage the community in working toward an ever improved ElmerFEM
- Enable a sustainable system for ElmerFEM maintainers to review contributions

Please follow these steps to have your contribution considered by the maintainers:

1. Always make sure your forked branch is derived from a more or less up-to-date `devel` branch
1. Follow the [styleguides](#styleguides)
3. Run all ´ctest´ to verify that the code enhancement do not affect the behaviour of other parts of the code


While the prerequisites above must be satisfied prior to having your pull request reviewed, the reviewer(s) may ask you to complete additional design work, tests, or other changes before your pull request can be ultimately accepted.

## Styleguides

### F90 style guide

Elmer uses 2008 Fortran standard in which all new code should be provided.

for formatting we use the following rules

* we use 2 space indent
* names functions/subroutines are written in Pascal case (or Upper Camel Case), e.g. OptimizeBandwidth()
* variable names for short variables are small caps, while compound names usually in Pascal case also but there some variation. 
* Fortran-language specific elements, such as variable-types (e.g., INTEGER) and properties (e.g., ALLOCATABLE) are written in capital letters
* Comments inside the code 

### Git Commit Messages

* Use the present tense ("Add feature" not "Added feature")
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
* Limit the first line to 72 characters or less
* Reference issues and pull requests liberally after the first line
* When only changing documentation, include `[ci skip]` in the commit title

### Documentation Styleguide

ElmerFEM model and solver documentation is written in LaTeX and maintained by CSC outside this repository. Part of the documentation is publicly
available in https://github.com/ElmerCSC/elmerfem-manuals. There are certain sub-projects (like Elmer/Ice) that provide their own documentation within this repository and use markdown language.


## Additional Notes

### Issue and Pull Request Labels

This section lists the labels we use to help us track and manage issues and pull requests. 

[GitHub search](https://help.github.com/articles/searching-issues/) makes it easy to use labels for finding groups of issues or pull requests you are interested in. We encourage you to read about [other search filters](https://help.github.com/articles/searching-issues/) which will help you write more focused queries.

The labels are loosely grouped by their purpose, but it's not required that every issue has a label from every group or that an issue can't have more than one label from the same group.

Please open an issue on `elmerfem/devel` if you have suggestions for new labels, and if you notice some labels are missing on some repositories, then please open an issue on that repository.

#### Type of Issue and Issue State


| Label name | `elmerfem/devel` :mag_right:| Description |
| --- | --- | --- |
| `enhancement` | [search][search-elmerfem-repo-label-enhancement]  | Feature requests. |
| `bug` | [search][search-elmerfem-repo-label-bug] | Confirmed bugs or reports that are very likely to be bugs. |
| `question` | [search][search-eerlmfem-repo-label-question] |  Questions more than bug reports or feature requests (e.g. how do I do X). |


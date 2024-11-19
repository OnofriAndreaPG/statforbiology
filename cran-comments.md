# Version 0.9.9

* This is the re-submission of the same statforbiology package

## Fixing the following comments from reviewer and resubmitting

* From: Konstanze Lauseker <konstanze.lauseker@wu.ac.at>
* Date: 02.07.2024 17:48

  - Please rather use the Authors@R field and declare Maintainer, Authors and Contributors with their appropriate roles with person() calls.
    * Changed. Thank you.
  - Please do not start the description with "It", "This package", package name, title or similar.
    * Changed. Thank you
  - Please write references in the description of the DESCRIPTION file in the form authors (year) <doi:...> .. -> please add more information about the blog.
    * Changed. Thank you
  - Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar)
    * Corrected. Thank you
  - You still have the standardtexts in your .Rd-files. Please omit them since they are not necessary.
    * Removed. Thank you
  - You write information messages to the console that cannot be easily suppressed. It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object. Instead of print() rather use message()/warning() or if(verbose)cat(..) (or maybe stop()) if you really have to write text to the console. (except for print, summary, interactive functions)
    * Corrected. Thank you
  - Please do not modify the user's global environment or the user's home filespace in your examples or vignettes by deleting objects rm(list = ls()) -> man/anova.aovList.Rd, man/SSexpoDecay.Rd
    * Command removed. Thank you
  
## Test environments

* local R installation, R 4.4.0
* win-builder (release)

## R CMD check results

0 errors | 0 warnings | 0 notes



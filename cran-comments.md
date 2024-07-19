## Resubmission

This resubmission incorporates the following revisions based on the reviewer's feedback:

* I have utilized single quotes for the package name in the description file.
* Corrected the citation formatting in the description file. 
* Replaced \dontrun with \donttest.
* I fixed all documentation issues.
* I made sure not to change user's options, par or working directory.

Additionally, significant enhancements have been made to the package documentation.

## R CMD check results

0 errors | 0 warnings | 2 note

* This is a new release.

❯ checking installed package size ... NOTE
    installed size is  6.0Mb
    sub-directories of 1Mb or more:
      extdata   5.0Mb
      
* The package includes raster files essential for demonstrating the 
  functionality of its features. These files are necessary as they are 
  extensively used throughout the documentation.
  
### Test environments

- R-hub windows-x86_64-devel (r-devel)
- R-hub  macOS arm64 (R-devel)
- R-hub ubuntu-clang
- R-hub  macOS (R-devel)
- R-hub  windows (R-devel)
- windows-latest (release; on GitHub Actions)
- macOS-latest (release; on GitHub Actions)
- ubuntu-latest (release; on GitHub Actions)
- ubuntu-latest (devel; on GitHub Actions)
- ubuntu-latest (older-1; on GitHub Actions)
- local ubuntu 22.04.4, install, R 4.4.1
- local OS X 14.2 install, R 4.3.2
- local windows 10 pro install, R 4.3.1

Function demo for loopBouquetPlotDocumentation
================
Authors: Jianhong Ou[^1]<br/>
Last modified: 2024-02-07

## Overview

The `loopBouquetPlot` is a new method introduced by **trackViewer**
package to visualize genomic interactions along with annotation for NGS
dataset such as HiC, HiChIP, PLAC-seq, ChIA-PET, and HiCAR data.

### Pre-requisites

- Basic knowledge of R syntax
- Basic knowledge of Docker
- Basic knowledge of shell commands
- A computer with internet connection

### Installation

To install this package, start R and enter:

``` r
library(BiocManager)
BiocManager::install("jianhong/loopBouquetPlotDocumentation")
```

### Documentation

To view documentation of loopBouquetPlotDocumentation, start R and
enter:

``` r
browseVignettes("loopBouquetPlotDocumentation")
```

### Contributions and Support

If you would like to contribute to this package, the standard workflow
is as follows:

1.  Check that there isn’t already an issue about your idea in the
    [jianhong/loopBouquetPlotDocumentation/issues](https://github.com/jianhong/loopBouquetPlotDocumentation/issues)
    to avoid duplicating work. If there isn’t one already, please create
    one so that others know you’re working on this
2.  [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo)
    the
    [jianhong/loopBouquetPlotDocumentation](https://github.com/jianhong/loopBouquetPlotDocumentation)
    to your GitHub account
3.  Make the necessary changes / additions within your forked repository
    following [Bioconductor
    contribution](https://contributions.bioconductor.org/)
4.  Use `devtools::build` and `devtools::check` to check the package
    work properly.
5.  Submit a Pull Request against the `master` or current
    `RELEASE_VERSION` branch and wait for the code to be reviewed and
    merged

If you’re not used to this workflow with git, you can start with some
[docs from
GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests)
or even their [excellent `git` resources](https://try.github.io/).

For further information or help, don’t hesitate to get in touch on the
[Bioconductor support site](https://support.bioconductor.org/) with tag
`#trackViewer`.

### Reporting bug/issues

Many thanks for taking an interest in improving this package. Please
report bug/issues at
[jianhong/loopBouquetPlotDocumentation/issues](https://github.com/jianhong/loopBouquetPlotDocumentation/issues).

[^1]: Regeneration Center, Duke University, Durham, North Carolina, USA.

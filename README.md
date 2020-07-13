## Instructions for getting started.


### Get the code

    cd ~
	git clone git@github.com:sje30/chickcones.git
	cd chickcones
	

### Install the necessary packages

Start R and then install the following:

	R
	install.packages(c("readxl", "devtools"))
	devtools::install_github("sjedrp")
    devtools::install_github("sjevor")
    devtools::install_github("sjedist")
    devtools::install_github("sjedmin")
    devtools::install_github("sjedrp")


### Run the code

First edit the script to run a field, e.g. field = "DT1" with the editor:

    gedit runone.R
	
(Or use Rstudio to change the field variable in that file.)

And then from the R prompt:

    source('runone.R')
	



	
    
	

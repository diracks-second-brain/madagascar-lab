# Makefile
# 
# Purpose: Automate building.
# 
# Site: https://www.geofisicando.com
# 
# Version 1.0
# 
# Programmer: Rodolfo A C Neves (Dirack) 05/10/2020
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

MAIN_SRC := model2d_smooth
DEP_SRC := velocity

MAIN := $(addsuffix .c, $(MAIN_SRC))
MAIN_EXEC := $(addsuffix .x, $(MAIN_SRC))
DEP := $(addsuffix .h, $(DEP_SRC))

RSF_INCLUDE := -I$$RSFSRC/include
RSF_LIB := -L$$RSFSRC/lib
RSF_LINK := -lrsf -lm -lblas

all:	$(MAIN_EXEC)

%.x:	%.c $(DEP)
	gcc $< $(RSF_INCLUDE) $(RSF_LIB) -o $@ $(RSF_LINK)

clean:
	rm *.x
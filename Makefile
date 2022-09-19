INSTALL_DIR = ~/unit-apps/ClassPro/0.5/

SRC = src
BIN = bin
SCRIPTS = scripts

ALL = ClassPro ClassGS prof2class class2acc class2cns

all:
	mkdir -p $(BIN) && cd $(SRC) && make && cp $(ALL) ../$(BIN)

clean:
	rm -f $(ALL) && cd $(SRC) && make clean

install:
	mkdir -p $(INSTALL_DIR)
	cp -r $(BIN) $(SCRIPTS) $(INSTALL_DIR)

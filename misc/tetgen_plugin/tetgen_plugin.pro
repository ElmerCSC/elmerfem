ELMER_HOME=$$(ELMER_HOME)

isEmpty(ELMER_HOME) {
   error("ELMER_HOME is undefined")
}

message(ELMER_HOME=$${ELMER_HOME})

TEMPLATE = subdirs
SUBDIRS = plugin testapp

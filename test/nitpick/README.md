# Nitpick

## Files

* `exprcase.hh` - Defines the ExprCase class used in the test cases. It is referenced in every test/case/*_def.hh file
* `nitpick_include.hh` - Common includes and macros, shared by both .cc files. Throws preprocessor error if DEF_FILE or GEN_FILE is not defined
* `nitpick_generate.cc` - Takes the provided DEF_FILE, uses DagPrinter to generate the output and writes to GEN_FILE
* `nitpick_run.cc` - Takes the provided DEF_FILE and GEN_FILE and runs and times the expression for different vec_sizes
MSASolver
=========

Heuristic search for Multiple Sequence Alignment (MSA)

This is project contains an implementation of the MSA domain with
affine gap costs and a version of the A* algorithm for computing
minimum cost alignments.  It's provided as a Maven project written in
Java.  Below are instructions for building and running a simple
experiment.

### Build

1. Install a recent version of Maven
2. Run 'mvn install' in the same directory as pom.xml

### Run

1. Install a recent version of Ant
2. Run a simple experiment with 'ant -f msa.xml'

### Experiments

When you run as described above the solver runs with the A* algorithm
and solves a very simple instance.  Statistics are output to a
plain text file.  The msa.xml defaults to this simple setup.  To run
the solver on the full set of "easy" instances (solvable by A* on a
machine with at least 10GB) run the following command:

    'ant -f msa.xml MSASolver.balibase'

You can specify an alternative output directory with -Doutput

    'ant -f msa.xml MSASolver.balibase -Doutput="/my/output/directory"'

You can specify an alternative input directory with -Dinput

    'ant -f msa.xml MSASolver.balibase -Doutput="/my/input/directory"'

Input files must end in '.seq' and be formatted the same way the
included seq files are formatted.

### Files

msa/ref1_seq - a directory containing files for BAliBASE Ref.1 benchmarks
msa/ref1_seq_easy - a directory containing files for BAliBASE Ref.1 benchmarks
msa/pam250.sub - a file that contains the PAM250 cost matrix
src - a  directory that contains the Java src files
msa.xml - an ant script for running experiments
pom.xml - Maven POM file for building and collecting dependencies
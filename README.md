
This is the companion code to the paper https://hal.inria.fr/hal-01522903/

To compile, you need to install cgal and sqlite3.

Then mkdir build; cd build; cmake ..; make

Then you can play with the executable `tettest`.

Or start the benchmark with the Python script.

Data will be written to a  SQLite database.

The other Python script can read this data and output the plots, as seen in the paper.

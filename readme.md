# Compile
Using cmake files to compile the source files. For instance,
```
cd biclique
cmake .
make
```

# Usage
To execute the code, you need to run the following executable files, which accept the following optional parameters:

- "-l": The left size constraint.

- "-r": The right size constraint.

- "-k": the value of k for k-plex, 

- "-d": Selceted from 'two' or 'core'.

An running example:

```
./biclique/bin/MBC [filepath] -a=[1,2] -k=[0,n] -e=[0,1]
```

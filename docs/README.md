# Documentation
This folder will contain any documentation that might make it simpler to
parse and contribute to this project. We will maintain the overall architecture
here and simple high-level descriptions of the key data strucures. 

## The G-Tree
The "G-Tree" is simply a trie data structure containing a series of motifs,
indexed on valid base pair transitions. The "central conceit" of our 
achitecture is that clever allocation and management of resources can 
obviate most of the need for complex data structures and algorithms currently 
central to read alignment.


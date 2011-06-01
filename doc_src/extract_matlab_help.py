#!/usr/bin/env python
import re
import os

def matlab2asciidoc(matlab_file, path='./funcs'):

    # Grab the text from the first comment block
    helptext = []
    gothelp = False
    inhelp = False
    fmat = open(matlab_file, 'r')
    for line in fmat:
        if gothelp:
            exit
        elif re.match(r'\s*\%', line):
            inhelp = True
            helptext.append(re.sub(r'\s*\%', '', line))
        elif inhelp:
            gothelp = True

    fmat.close()

    # Open an output file:
    adoc_name = os.path.basename(os.path.splitext(matlab_file)[0])
    adoc_name = adoc_name + '.txt'
    fadoc = open(os.path.join(path, adoc_name), 'w')

    #Function name must match file name
    fadoc.write('== ' +
        os.path.basename(os.path.splitext(matlab_file)[0]) + '\n\n')
    for line in helptext:
        fadoc.write(line) 
    fadoc.close()

def asciidoc_index(files, path='./funcs'):

    index_name = 'index.txt'
    find = open(os.path.join(path, index_name), 'w')
    find.write('== MSAT function documentation listi\n\n')
    for file in files:
        fname = os.path.basename(os.path.splitext(file)[0])
        find.write('* link:./' + fname + '.html[' + fname + ']\n')

    find.close()
        
    

if __name__ == '__main__':
    import sys
    for filename in sys.argv[1:]:
        matlab2asciidoc(filename)

    asciidoc_index(sys.argv[1:])

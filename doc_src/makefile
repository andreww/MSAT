default: html_funcs html_docs

AD=asciidoc 

asciidoc_funcs:
	(mkdir -p funcs; ./extract_matlab_help.py ../msat/*.m)

html_funcs: asciidoc_funcs
	(mkdir -p funcs; cd funcs ; for f in *.txt ; do $(AD) $$f ; done)

html_docs: MSATstart.txt MSATguide.txt MSATfunctions.txt ../RELEASE_NOTES
	$(AD) -a latexmath MSATstart.txt
	$(AD) MSATguide.txt
	$(AD) -a latexmath MSATfunctions.txt
	$(AD) -o MSATrelease_notes.html ../RELEASE_NOTES 
	$(AD) -o MSATreadme.html ../README

install: html_funcs html_docs
	mkdir -p ../msat/html/funcs
	cat helptoc.xml_header helptoc.xml_funcs helptoc.xml_footer > ../msat/html/helptoc.xml
	cp MSATrelease_notes.html MSATreadme.html MSATfunctions.html  MSATstart.html MSATguide.html albite.png MSATfunctions.html ../msat/html
	cp ./funcs/*.html ../msat/html/funcs
	$(MAKE) clean

external-install: html_funcs html_docs
	# NB - $(TARGET) must be set on the command line 
	# or this will eat your computer.
	mkdir -p $(TARGET)/funcs
	cp MSATstart.html MSATguide.html MSATreadme.html MSATrelease_notes.html MSATfunctions.html albite.png $(TARGET)
	cp ./funcs/*.html $(TARGET)/funcs
	$(MAKE) clean

clean:
	rm -f ./funcs/*.html ./funcs/*.txt helptoc.xml_funcs
	rm -f MSATstart.html MSATguide.html MSATfunctions.html
	rm -f MSATrelease_notes.html MSATreadme.html
	rm -rf funcs

clean_all: clean
	rm -rf ../msat/html/*	


LDDist-1.3.2
============

LDDist-1.3.2 - Hopefully running on Ubuntu. All credit to original authors.

Original content from http://artedi.ebc.uu.se/molev/software/LDDist.html

I'm putting this here as I think it is a better way of controlling changes - even minor bug fixes - rather than just adding it to a long forgotten blog post.

Compiling generally (with make -v = 3.81 and gcc -v = 4.6.3 on Ubuntu 10.4+) gives the below errors

LDDist_wrap.cxx:442:9: error: expected unqualified-id before string constant
LDDist_wrap.cxx:443:9: error: ‘SwigPerlWrapper’ does not name a type
LDDist_wrap.cxx:448:3: error: ‘SwigPerlWrapperPtr’ does not name a type
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx:1167:1: error: too many initialisers for ‘swig_command_info’
LDDist_wrap.cxx: In function ‘void boot_LDDist(CV*)’:
LDDist_wrap.cxx:1186:62: error: ‘struct swig_command_info’ has no member named ‘wrapper’
make: *** [LDDist_wrap.o] Error 1


But changing line 442 from:

typedef XS(SwigPerlWrapper);

to

typedef XSPROTO(SwigPerlWrapper);

seems to allow it to compile properly.

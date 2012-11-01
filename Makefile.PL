#File : Makefile.PL
use ExtUtils::MakeMaker;
WriteMakefile (
	'NAME'	=>	'LDDist',
	'LIBS'	=>	['-lm -lstdc++'],
	'OBJECT'	=>	'LDDist.o LDDist_wrap.o'
);
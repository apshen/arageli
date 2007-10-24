@rem *****************************************************************************
@rem
@rem     create_guide.bat
@rem
@rem     This file is a part of the Arageli library.
@rem
@rem     Copyright (C) 2005--2007 Sergey S. Lylain
@rem     University of Nizhni Novgorod, Russia
@rem
@rem     The Arageli Library is free software; you can redistribute it and/or
@rem     modify it under the terms of the GNU General Public License version 2
@rem     as published by the Free Software Foundation.
@rem
@rem     The Arageli Library is distributed in the hope that it will be useful,
@rem     but WITHOUT ANY WARRANTY; without even the implied warranty of
@rem     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
@rem     GNU General Public License for more details.
@rem
@rem     You should have received a copy of the GNU General Public License along
@rem     with this program; if not, write to the Free Software Foundation, Inc.,
@rem     51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
@rem
@rem     We are also open for dual licensing for the whole library or
@rem     for its particular part. If you are interested to get the library
@rem     in this way, i.e. not under the GNU General Public License,
@rem     please contact Arageli Support Service support.arageli@gmail.com.
@rem
@rem *****************************************************************************

@rem Building "Arageli: a Library for Algebraic Computations" paper.

pdflatex ..\doc\src\arageli_overview_announce_preprint.tex
pdflatex ..\doc\src\arageli_overview_announce_preprint.tex
del arageli_overview_announce_preprint.aux
del arageli_overview_announce_preprint.log
move arageli_overview_announce_preprint.pdf ..\doc\

@rem Building "Arageli Coding Guidelines".

cd ..\doc\src\acg\rus
latex acg.tex
latex acg.tex
dvips acg.dvi
ps2pdf acg.ps
del acg.aux
del acg.log
del acg.toc
del acg.dvi
del acg.ps
del acg.out
move acg.pdf ..\..\..
cd ..\..\..\..\build

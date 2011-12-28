#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV, Andrey Simurzin
##


#----------------------------------------------------------------------------
from Foam import ref, man

#----------------------------------------------------------------------------
def _createFields( runTime, mesh ):
    ref.ext_Info() << "Reading field T\n" << ref.nl
    T = man.volScalarField( man.IOobject( ref.word( "T" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                             mesh )

    ref.ext_Info() << "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )


    ref.ext_Info() << "Reading transportProperties\n" << ref.nl
    transportProperties = man.IOdictionary( man.IOobject( ref.word( "transportProperties" ),
                                                          ref.fileName( runTime.constant() ),
                                                          mesh,
                                                          ref.IOobject.MUST_READ_IF_MODIFIED,
                                                          ref.IOobject.NO_WRITE ) )


    ref.ext_Info() << "Reading diffusivity D\n" << ref.nl

    DT = ref.dimensionedScalar( transportProperties.lookup( ref.word( "DT" ) ) )

    phi = man.createPhi( runTime, mesh, U )
           
    return T, U, transportProperties, DT, phi 


#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )
    
    T, U, transportProperties, DT, phi = _createFields( runTime, mesh )

    simple = man.simpleControl( mesh )

    ref.ext_Info() << "\nCalculating scalar transport\n" << ref.nl
        
    CoNum, meanCoNum = ref.CourantNo( mesh, phi, runTime )
    
    while simple.loop():
        ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl
        
        while simple.correctNonOrthogonal():
            ref.solve( ref.fvm.ddt( T ) + ref.fvm.div( phi, T ) - ref.fvm.laplacian( DT, T ) )
            pass

        runTime.write()
        pass
        
    ref.ext_Info() << "End\n" << ref.nl 

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
import sys, os
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020100" ):
   if __name__ == "__main__" :
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   ref.ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.1.0 or higher \n "
   pass

	
#--------------------------------------------------------------------------------------


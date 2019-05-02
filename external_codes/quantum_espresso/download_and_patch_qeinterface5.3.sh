#!/bin/sh

# Attempt to automatically download and patch Quantum-Espresso for pw2qmcpack converter 
# Patch developed by William Parker / Argonne National Lab
# This simple script, Paul Kent / Oak Ridge National LaB

codename=espresso-5.3.0
archivename=qe-5.3.tar.gz
rootdirname=q-e-qe-5.3
if [ ! -e ${archivename} ]; then
echo --- Did not find ${archivename} in current directory.
# Full URL of espresso download link obtained from qe-forge on 29 September 2014
# Will need to be updated between versions
qeurl=https://github.com/QEF/q-e/archive/${archivename}
echo --- Attempting to download. This can take several minutes. 
wget  ${qeurl}
else
echo --- Found and using ${archivename} in current directory.
fi
if [ ! -e ${archivename} ]; then
echo --- ERROR: Could not find ${archivename}
echo --- Something went wrong... possibly a bad URL for the file download or you are offline
echo --- Please advise QMCPACK Developers via Google Groups if this problem persists
else
echo --- Unpacking
tar xvzf ${archivename}
if [ ! -e ${rootdirname}/PW/src/Makefile ]; then
echo --- ERROR: Could not find PW/src/Makefile
echo --- Something went wrong... probably a failure to download the full archive.
echo --- Check ${archivename}. Delete if a partial download and retry.
echo --- Also check $qeurl is valid - perhaps the files have moved online.
echo --- Please advise QMCPACK Developers via Google Groups if this problem persists
else
cd ${rootdirname}
patch -p1 -i ../add_pw2qmcpack_to_int${codename}.diff
cd ..
if [ -e $rootdirname/PW/src/pwinterfacemod.f90 ]; then
echo --- SUCCESS: int${codename} patched for pwscf C++ interface
echo "--- Configure using ./configure --with-hdf5 HDF5_DIR=(HDF5 base directory)"
echo --- Add platform specific options as needed
else
echo --- ERROR: Could not find PW/src/pwinterfacemod.f90 after patching
echo --- Probably the patch is missing or the archive has been updated.
echo --- Please advise QMCPACK Developers via Google Groups.
fi    
fi

fi


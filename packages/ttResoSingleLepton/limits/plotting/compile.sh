echo "Compiling TTBarSystematics.cpp..."
g++ `$ROOTSYS/bin/root-config --cflags --glibs` -o runSystematics TTBarSystematics.C
echo "...done"
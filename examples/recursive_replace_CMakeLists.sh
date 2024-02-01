#!/bin/bash  
  
# this bash script recursively replace all CMakeLists.txt in example subfolders with batch_CMakeLists.txt  
  
for f in $(find . -name 'CMakeLists.txt'); do  
    if [ -f "$f" ]; then  
        cp batch_CMakeLists.txt "$f"  
        echo replaced "$f"  
    else  
        # file doesn't exist  
        echo "File does not exist: $f"  
    fi  
done  

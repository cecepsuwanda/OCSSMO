#!/bin/bash

./ocssmo ~/Dataset/NSL/model_100_2_svm_29attr_100 ~/cpp/my_svm/data/kddcupSVM_32attr.names 0.0001 0.0001 0.5 0.5 0.9 0.9 0.01 0.01 | tee hasil/hasil.txt

#!/bin/sh
echo -n "Is conda installed? : (y/n)"
read conda_installed

if [ $conda_installed = "n" ]; then
   echo 
else
    wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    bash Miniconda3-latest-Linux-x86_64.sh
    rm Miniconda3-latest-Linux-x86_64.sh


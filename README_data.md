# README_data.md

## Introduction

Documentation for FreeDV data channel.
A work in progres..

## Quickstart

1. Simple test using mode 2400A

  ```
   $ ./src/freedv_data_tx 2400A - --frames 15 | ./src/freedv_data_rx 2400A -  

  ```

2. Same for 2400B and 800XA

  ```
   $ ./src/freedv_data_tx 2400B - --frames 15 | ./src/freedv_data_rx 2400B -  

  ```

  ```
   $ ./src/freedv_data_tx 800XA - --frames 15 | ./src/freedv_data_rx 800XA -  

  ```

3. Using a different callsign and secondary station id

  ```
   $ ./src/freedv_data_tx 2400A - --callsign T3ST --ssid 15 --frames 15 | src/freedv_data_rx 2400A -  
  ```

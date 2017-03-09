Bandeirantes: A Graph-based Approach for Curve Tracing and Boundary Tracking
===================

We provide the supplementary materials for the following scientific publication of the Bandeirantes method, which extends the previous method live wire on the fly for curve tracing and boundary tracking.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project.

### Prerequisites

* Libraries:
    * [munkres-cpp](https://github.com/saebyn/munkres-cpp) - Hungarian algorithm 
    * gft (in code)
        - optional: zlib: zlib1g-dev

* Images:
    * [Images for test the methods](http://www.vision.ime.usp.br/~mtejadac/bandeirantes.html) - Scanned completed answers sheets used in the experiment


## Build

### GFT
```bash
$ cd lib/gft
$ make clean && make
```
### Bandeirantes method
```bash
$ export GFT_DIR=$(pwd)/lib/gft
$ make clean && make
```

## Configurations files

### Format of configuration of each image
<pre>
n
x<sub>1</sub>    x<sub>1</sub>     [solution]
x<sub>2</sub>    x<sub>2</sub>     [solution]
...
x<sub>n</sub>    x<sub>n</sub>     [solution]
m
x<sub>1</sub>    x<sub>1</sub>   
x<sub>2</sub>    x<sub>2</sub>   
...
x<sub>m</sub>    x<sub>m</sub> 
</pre>
where: 
<pre>n: number of initial points
m: number of terminal points
[solution]: solution corresponding to order in list m 
</pre>


Example:
```
7
148   99  7
139  165  4
129  205  6
123  246  2
111  315  5
103  379  3
 90  450  1
7
483   88
466  160
441  228
427  291
408  352
390  403
372  462
```

![image_test.pgm](http://www.vision.ime.usp.br/~mtejadac/bandeirantes_files/graph0002a.png)

### Configuration file of experiments table
Image list file
<pre>/path_to_directory_images/
image_001.pgm
image_002.pgm
image_003.pgm
...
image_last.pgm
</pre>
 
 Configuration file
<pre>/path_to_configuration_files/
conf_image_001.txt
conf_image_002.txt
conf_image_003.txt
...
conf_image_last.txt</pre>

## Run

### Run individual image
Usage:
```
main <image_file> <config_file> <method>
method: 0 ... Riverbed
        1 ... LiveWire
        2 ... G-wire
        3 ... Bandeirantes
```

Example:
```bash
$ ./main database/img_001.pgm database/img_001.txt 3
```

### Run list images
Usage:
```
runexp <image_list> <config_list> <method>
method: 0 ... Riverbed
        1 ... LiveWire
        2 ... G-wire
        3 ... Bandeirantes
```

Example:
```bash
$ ./runexp txt/imglist.txt txt/conflist.txt 3
```
## Authors

* **Marcos A. Tejada Condori** - [website](http://www.vision.ime.usp.br/~mtejadac/)

* **Paulo A. V. de Miranda** - [website](http://www.vision.ime.usp.br/~pmiranda/)

* **Lucy Choque Mansilla** - [website](http://www.vision.ime.usp.br/~lucyacm/)

## License

This project is licensed under the terms of the MIT License.
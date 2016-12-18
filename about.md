# Numerical Integration (CNI)

Newton-Cotes formulas are not accurate for a large integral domain. The main idea of the CNI is to divide integral domain into subdomains and use simple integration rule for each subdomain. 

Using the C++ code (composite.cpp), I computed composite numerical integration over a quadrilateral region. We described the composite numerical integration over an arbitrary quadrilateral region in our following paper: 

Md. Shafiqul Islam and M. Alamgir Hossain (2009)- Numerical Integrations over an Arbitrary Quadrilateral Region, Applied Mathematics and Computations Volume 210, Issue 2, 15 April 2009, Pages 515-524.

## Run
Example 2 of Islam and Hossain (2009) can run in a linux machine with the following instructions: 

$ g++ composite.cpp<br>
$ ./a.out

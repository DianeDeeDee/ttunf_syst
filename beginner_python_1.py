import random
import sys
import os


def factorial(n):

    """factorial() is the function that will calculate the
    the factorial of a non-negative integer n
    """
    res = 1
    for i in range(1, n+1):
        res=res*i
    return res

#print(factorial(3),"=3!")
#print(factorial(5),"=5!")
#print("\n\n")

def is_prime(number):
    """is_prime() is the the function that tests the primality of the input number will check if number is prime.
    If it's prime, boolean is_it_prime is True, otherwise it' is False.
    In number theory, Wilson's theorem states that a natural number > 1 is a prime number
    if and only if ( number − 1 ) !   ≡   − 1 ( mod number )
    Let's define:
      fact = ( number − 1 ) ! + 1
      mod = fact % number
        """

    fact= (factorial(number - 1) + 1)
    mod = fact % number
    for number in range(2,1000000):
        if (mod == 0):
            is_it_prime = True
        else:
            is_it_prime = False

    return is_it_prime

#print("factorial(b-1)=",factorial(-1),"; factorial(b-1)+1=",factorial(5-1)+1," "
 #      ";fac%b=",(factorial(5-1)+1)%5)
print("\n\n")
print("- Come on Magic Ball give me the answer! Is %s a prime? "%"29")
print("-",is_prime(29),"!")
print("- How cool it's working! Let's try another number...")

print("\n\n")

print("We can improve the code so it will check all numbers in range, and not "
      "just one specific number.")
print("Then, we can using an output.txt where all prime number are saved automatically. ;)")

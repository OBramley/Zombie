import hello
print(hello.__doc__)

print(hello.logical_to_integer.__doc__)

sieve_array=hello.sieve(100)
prime_numbers=hello.primes.logical_to_integer(sieve_array, sum(sieve_array))
print(prime_numbers)

    unsigned long int block_size = 1048576;
    // unsigned long int block_size = 2;
    unsigned long long int block_low_value;
    unsigned long long int block_high_value;
    for (block_low_value = low_value; block_low_value <= high_value; block_low_value += 2 * block_size) {
        block_high_value = block_low_value + 2 * (block_size - 1);
        if (block_high_value > high_value)
            block_high_value = high_value;

        index = 0;
        prime = 3;
        while (prime * prime <= block_high_value) {
            if (prime * prime > block_low_value)
                first = (prime * prime - block_low_value) / 2;
            else {
                if (!(block_low_value % prime))
                    first = 0;
                else if (block_low_value % prime % 2 == 0)
                    first = prime - ((block_low_value % prime) / 2);
                else
                    first = (prime - (block_low_value % prime)) / 2;

            }
            first += (block_low_value - low_value) / 2;
            for (i = first; i < (block_high_value - low_value) / 2; i += prime) {
                marked[i] = 1;
            }
            while (local_prime_marked[++index]);
            prime = index * 2 + 3;
        }
    }

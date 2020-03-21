NAME
    cudaKeccakSolver

DESCRIPTION
    C++/CUDA implementation of finding the 4-round Keccak preimage on parameter [r=40, c=160, nr=4]

DOCUMENTATION
    implementation detail is on preparation

DEPENDENCIES
    cmake, make, gcc, CUDA

BUILD
    $ mkdir build && cd build
    $ cmake .. && make

USAGE
    $ ./keccaksolver --t [cpu_thread_num] --gb_start [guessing_bits_start] --gb_end [guessing_bits_end] --dev_id [gpu_device_id]

EXAMPLE
    $ ./keccaksolver --t 20 --gb_start 0x12345 --gb_end 0xabcde --dev_id 0

        start Keccak solver with 20 cpu threads, guessing bits starts at 0x12345 and ends at 0xabcde
        using gpu device 0 to solve mq sub systems

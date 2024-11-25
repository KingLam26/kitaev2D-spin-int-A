"""
This is the main Python file to execute. It utilizes parallel processing to handle large-scale computations for spin systems involving permutations
and energy summations. It employs producer-consumer architecture for efficient processing and uses multiprocessing 
to distribute the workload across multiple CPU cores.

Key Features:
1. **Parallel Processing**:
    - A producer generates all unique permutations of a given bond configuration and places them in a shared queue.
    - Multiple consumers process the permutations in parallel, calculating energy or coefficients and updating shared results.

2. **Shared Memory Management**:
    - Results and counters are stored in shared memory (`mp.Value`) to enable safe concurrent access.
    - Locks (`mp.Lock`) are used to manage access to shared variables, ensuring thread safety.

3. **Dynamic Runtime Feedback**:
    - Periodically prints progress updates, including the percentage of permutations processed, runtime, estimated time remaining, and processing speed.

4. **Graceful Termination**:
    - Handles `KeyboardInterrupt` gracefully, ensuring processes terminate cleanly.
    - Sends sentinel values (`None`) to signal consumers when no more tasks are available.

5. **Modifications and Scalability**:
    - Scales across multiple cores using a configurable number of consumer processes.
    - Allows runtime modifications to summation factors (e.g., spin factors, multiplicative adjustments from research papers).

6. **Main Execution Flow**:
    - `main()` orchestrates the initialization, producer-consumer setup, and final result aggregation.
    - Input arguments define the bond permutation configuration and number of cores to utilize.

7. **Runtime Statistics**:
    - Tracks start and end times, total runtime, and provides detailed performance metrics.
    - Results, including runtime and computed sums, are saved for further analysis.
"""

from core import *
from utils import *
from datetime import datetime
import multiprocessing as mp
import sys, time


##### parallel processing functions #####
def consumer_task(queue, results_values, value_locks, counter, counter_lock, tot_bond_permutes, start_time):
    while True:
        bond_permute = queue.get()
        if bond_permute is None:
            # Sentinel received; exit consumer
            break

        # Process the permutation
        if compute_switch:
            result = process_permute(list(bond_permute))[1]
        #print(bond_permute, result)

        # Safely add the result to one of the shared lists
        list_index = counter.value % len(results_values)  # Round-robin assignment to lists
        with value_locks[list_index]:
            current_sum = results_values[list_index].value
            new_sum = current_sum + result
            results_values[list_index].value = new_sum

        # Safely update the counter and check if it's a multiple of 1000
        with counter_lock:
            counter.value += 1
            if counter.value % num_permutes_print == 0:
                now_time = datetime.now()
                runtime = (now_time - start_time).total_seconds() / 60
                num_permutes_left = tot_bond_permutes - counter.value
                time_left = (runtime / counter.value) * num_permutes_left
                percent = round((counter.value / tot_bond_permutes) * 100, 2)
                its_per_sec = counter.value / (runtime * 60)
                queue_size = queue.qsize()
                print(f"{counter.value}<{tot_bond_permutes}, {percent}%, {round(runtime, 1)}<{round(time_left, 2)} mins, {int(its_per_sec)} its/s, qsize: {queue_size}.")

def producer_task(queue, bond_permute, cores):
    # Generate permutations and put them in the queue
    for perm in gen_next_permute(bond_permute):
        queue.put(tuple(perm))

    # Send the sentinel (None) to stop the consumers
    for _ in range(cores):
        queue.put(None)


##### main script call main #####
def main():
    ##### process run script input params #####
    if len(sys.argv) > 2:
        print(f"Received input: {sys.argv[1]}, use {sys.argv[2]} cores")
        arg = sys.argv[1]
        cores = int(sys.argv[2])
        bp_permute_run = bp_dict[arg]
        tot_bond_permutes = num_bond_permutes(bp_permute_run)
        print(f"bp_permute_run: {bp_permute_run}, {tot_bond_permutes} permutations")

    else:
        print("insufficient input provided")    

    
    start_time = datetime.now()
    ##### parallel processing code block #####
    # Shared queue, shared total sum, shared counter
    queue = mp.Queue(maxsize=queue_maxsize)

    # Create several lockable lists
    results_values = [mp.Value('d', 0.0) for _ in range(num_lockable_lists)]
    value_locks = [mp.Lock() for _ in range(num_lockable_lists)]
    counter = mp.Value('i', 0)    # Shared counter for tracking progress, initialized to 0

    # Single lockable counter
    counter_lock = mp.Lock()

    # Define producers and consumers
    producer = mp.Process(target=producer_task, args=(queue, bp_permute_run, cores))
    consumers = [mp.Process(target=consumer_task, args=(queue, results_values, value_locks, counter, counter_lock, tot_bond_permutes, start_time)) for _ in range(cores)]
    
    try:
        # Start producer
        producer.start()

        # Start consumers
        for c in consumers:
            c.start()

        # Join producer
        producer.join()

        # Join consumers
        for c in consumers:
            c.join()

    except KeyboardInterrupt:
        print("KeyboardInterrupt detected, terminating processes...")

        producer.terminate()
        for c in consumers:
            c.terminate()

        time.sleep(1)

        producer.join()
        for c in consumers:
            c.join()

    ##### finish #####
    finally:
        # Summing Up Results Using Decimal
        print(f"Total permutations processed: {counter.value}, summing up results...")
        total_sum = sum(
            results_value.value  # Extract the value from each mp.Value
            for results_value in results_values
        )

        # modify total sum: spin factors of 2
        add_spin_factor = 2**(spin_S*8)
        total_sum/=add_spin_factor

        # modify total sum: any multiplicative factors from Koga's paper for easy comparison
        total_sum*=multi_factor

        print(f"Total sum of processed permutations: {total_sum}")

        # run time statistics
        end_time = datetime.now()
        runtime = end_time - start_time
        runtime_minutes = runtime.total_seconds() / 60

        start_time_str = start_time.strftime('%Y-%m-%d %H:%M:%S')
        end_time_str = end_time.strftime('%Y-%m-%d %H:%M:%S')
        
        # save results
        bp_label = arg
        save_final_results(bp_label, start_time_str, end_time_str, runtime_minutes, cores, total_sum)


##### execute #####
if __name__ == '__main__':
    # This guard is necessary for multiprocessing to work on Windows
    mp.freeze_support()
    main()
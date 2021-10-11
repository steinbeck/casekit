package casekit.threading;

import java.util.Collection;
import java.util.concurrent.*;
import java.util.function.Consumer;

public class MultiThreading {

    public static ExecutorService initExecuter(final int nThreads) {
        return Executors.newFixedThreadPool(nThreads);
    }

    public static void stopExecuter(final ExecutorService executor, final long seconds) {
        executor.shutdown();
        try {
            if (!executor.awaitTermination(seconds, TimeUnit.SECONDS)) {
                System.err.println("killing non-finished tasks!");
                executor.shutdownNow();
            }
        } catch (final InterruptedException e) {
            System.err.println("killing non-finished tasks!");
            executor.shutdownNow();
        }
    }

    public static <T> void processTasks(final Collection<Callable<T>> callables, final Consumer<T> consumer,
                                        final int nThreads, final long seconds) throws InterruptedException {
        // initialize an executor for parallelization
        final ExecutorService executor = initExecuter(nThreads);
        // execute all task in parallel
        executor.invokeAll(callables)
                .stream()
                .map(future -> {
                    try {
                        return future.get();
                    } catch (final InterruptedException | ExecutionException e) {
                        throw new IllegalStateException(e);
                    }
                })
                .forEach(consumer);
        // shut down the executor service
        stopExecuter(executor, seconds);
    }
}

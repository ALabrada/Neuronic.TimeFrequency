namespace Neuronic.TimeFrequency
{
    /// <summary>
    /// Abstraction of a stop criteria for an iterative algorithm.
    /// </summary>
    /// <typeparam name="T">The type of the algorithm's state.</typeparam>
    public interface IStopCriteria<in T>
    {
        /// <summary>
        /// Determines if the algorithm should stop iterating.
        /// </summary>
        /// <param name="state">The current state.</param>
        /// <returns><c>true</c> if the algorithm should stop; otherwise, <c>false</c>.</returns>
        bool ShouldStop(T state);
    }
}
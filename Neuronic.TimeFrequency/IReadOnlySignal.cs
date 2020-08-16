using System.Collections.Generic;

namespace Neuronic.TimeFrequency
{
    /// <summary>
    /// Read-only interface of a <see cref="Signal{T}"/>
    /// </summary>
    /// <typeparam name="T">The type of the samples.</typeparam>
    /// <seealso cref="System.Collections.Generic.IReadOnlyList{T}" />
    public interface IReadOnlySignal<out T> : IReadOnlyList<T>
    {
        /// <summary>
        /// Gets the time of the first sample.
        /// </summary>
        double Start { get; }
        /// <summary>
        /// Gets the sampling rate.
        /// </summary>
        double SamplingRate { get; }
        /// <summary>
        /// Gets the sampling period.
        /// </summary>
        double SamplingPeriod { get; }
        /// <summary>
        /// Gets the time of the last sample.
        /// </summary>
        double End { get; }
    }
}
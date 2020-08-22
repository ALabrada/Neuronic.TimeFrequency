﻿using System.Collections.Generic;

namespace Neuronic.TimeFrequency
{
    /// <summary>
    /// Abstraction of a Time Frequency Representation (TFR).
    /// </summary>
    public interface ITimeFrequencyRepresentation
    {
        /// <summary>
        /// Gets the sampling period in the time domain.
        /// </summary>
        double SamplingPeriod { get; }

        /// <summary>
        /// Gets the amount of samples in the time domain.
        /// </summary>
        int SampleCount { get; }

        /// <summary>
        /// Gets the TFR magnitude for the specified time and frequency.
        /// </summary>
        /// <param name="offset">The sample offset in the time domain.</param>
        /// <param name="frequency">The frequency.</param>
        /// <returns>The TFR magnitude.</returns>
        double this[int offset, double frequency] { get; }
    }

    /// <summary>
    /// Abstraction of a Bilinear or Quadratic Time Frequency Representation (QTFR).
    /// </summary>
    /// <seealso cref="Neuronic.TimeFrequency.ITimeFrequencyRepresentation" />
    public interface IBilinearTimeFrequencyRepresentation : ITimeFrequencyRepresentation
    {
        /// <summary>
        /// Gets the amount of samples in the frequency domain.
        /// </summary>
        int FrequencyCount { get; }

        /// <summary>
        /// Gets the computed frequencies.
        /// </summary>
        IEnumerable<double> Frequencies { get; }

        /// <summary>
        /// Gets the TFR magnitude for the specified time and frequency.
        /// </summary>
        /// <param name="offset">The sample offset in the time domain.</param>
        /// <param name="frequencyIndex">The frequency index in <see cref="Frequencies"/>.</param>
        /// <returns>The TFR magnitude.</returns>
        double this[int offset, int frequencyIndex] { get; }
    }
}
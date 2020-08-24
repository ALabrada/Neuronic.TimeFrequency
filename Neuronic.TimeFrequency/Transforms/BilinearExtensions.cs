using System;
using System.Collections.Generic;
using System.Linq;

namespace Neuronic.TimeFrequency.Transforms
{
    /// <summary>
    /// Contains extension methods for <see cref="IBilinearTimeFrequencyRepresentation"/>.
    /// </summary>
    public static class BilinearExtensions
    {
        private static int ArgMin<T>(this IEnumerable<T> items) where T : IComparable<T>
        {
            using (var enumerator = items.GetEnumerator())
            {
                if (!enumerator.MoveNext())
                    return -1;

                var min = enumerator.Current;
                var minIndex = 0;
                for (int index = 1; enumerator.MoveNext(); index++)
                {
                    if (min.CompareTo(enumerator.Current) > 0)
                    {
                        min = enumerator.Current;
                        minIndex = index;
                    }
                }

                return minIndex;
            }
        }

        /// <summary>
        /// Finds the index of the computed frequency that is closest to the specified value.
        /// </summary>
        /// <param name="transform">The transform.</param>
        /// <param name="frequency">The frequency.</param>
        /// <returns>The closest index to <paramref name="frequency"/>.</returns>
        public static int FindClosestIndexOfFrequency(this IBilinearTimeFrequencyRepresentation transform, double frequency)
        {
            return transform.Frequencies.Select(f => Math.Abs(frequency - f)).ArgMin();
        }

        /// <summary>
        /// Enumerates the amplitude values associated to the specified frequency, along the time axis.
        /// </summary>
        /// <param name="transform">The transform.</param>
        /// <param name="frequency">The frequency.</param>
        /// <returns>The values with frequency <paramref name="frequency"/>.</returns>
        public static IEnumerable<double> EnumerateValueOfFrequency(this IBilinearTimeFrequencyRepresentation transform, double frequency)
        {
            var freqIndex = FindClosestIndexOfFrequency(transform, frequency);
            if (freqIndex < 0 || freqIndex >= transform.FrequencyCount)
                return Enumerable.Empty<double>();
            return transform.EnumerateValuesOfFrequencyAt(freqIndex);
        }

        /// <summary>
        /// Enumerates the amplitude values associated to the frequency at the specified index, along the time axis.
        /// </summary>
        /// <param name="transform">The transform.</param>
        /// <param name="index">The frequency index.</param>
        /// <returns>The values with the frequency at <paramref name="index"/>.</returns>
        public static IEnumerable<double> EnumerateValuesOfFrequencyAt(this IBilinearTimeFrequencyRepresentation transform, int index)
        {
            return Enumerable.Range(0, transform.SampleCount).Select(i => transform[i, index]);
        }

        /// <summary>
        /// Enumerates the computed coordinates in the time domain.
        /// </summary>
        /// <param name="transform">The transform.</param>
        /// <returns>The coordinates in the time axis.</returns>
        public static IEnumerable<double> EnumerateTimes(this IBilinearTimeFrequencyRepresentation transform)
        {
            return Enumerable.Range(0, transform.SampleCount)
                .Select(i => transform.StartTime + transform.SamplingPeriod * i);
        }

        /// <summary>
        /// Finds the offset of the time domain sample that is closest to the specified value.
        /// </summary>
        /// <param name="transform">The transform.</param>
        /// <param name="time">The time.</param>
        /// <returns>The closest offset to <paramref name="time"/>.</returns>
        public static int FindClosestOffsetOfTime(this IBilinearTimeFrequencyRepresentation transform, double time)
        {
            return (int) Math.Round((time - transform.StartTime) / transform.SamplingPeriod);
        }

        /// <summary>
        /// Enumerates the amplitude values associated to the specified time, along the frequency axis.
        /// </summary>
        /// <param name="transform">The transform.</param>
        /// <param name="time">The time.</param>
        /// <returns>The values with time <paramref name="time"/>.</returns>
        public static IEnumerable<double> EnumerateValuesOfTime(this IBilinearTimeFrequencyRepresentation transform, double time)
        {
            var offset = transform.FindClosestOffsetOfTime(time);
            if (offset < 0 || offset >= transform.SampleCount)
                return Enumerable.Empty<double>();
            return transform.EnumerateValuesOfTimeAt(offset);
        }

        /// <summary>
        /// Enumerates the amplitude values associated to the time sample at the specified offset, along the frequency axis.
        /// </summary>
        /// <param name="transform">The transform.</param>
        /// <param name="offset">The offset.</param>
        /// <returns>The values with the offset <paramref name="offset"/> in the time domain.</returns>
        public static IEnumerable<double> EnumerateValuesOfTimeAt(this IBilinearTimeFrequencyRepresentation transform, int offset)
        {
            return Enumerable.Range(0, transform.FrequencyCount).Select(i => transform[offset, i]);
        }

        /// <summary>
        /// Enumerates the computed amplitude along both axes.
        /// </summary>
        /// <param name="transform">The transform.</param>
        /// <returns>The amplitude values.</returns>
        public static IEnumerable<double> EnumerateValues(this IBilinearTimeFrequencyRepresentation transform)
        {
            return from i in Enumerable.Range(0, transform.SampleCount)
                from j in Enumerable.Range(0, transform.FrequencyCount)
                select transform[i, j];
        }
    }
}
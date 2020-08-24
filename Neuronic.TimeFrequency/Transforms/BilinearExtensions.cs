using System;
using System.Collections.Generic;
using System.Linq;

namespace Neuronic.TimeFrequency.Transforms
{
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

        public static int FindClosestIndexOfFrequency(this IBilinearTimeFrequencyRepresentation transform, double frequency)
        {
            return transform.Frequencies.Select(f => Math.Abs(frequency - f)).ArgMin();
        }

        public static IEnumerable<double> EnumerateValueOfFrequency(this IBilinearTimeFrequencyRepresentation transform, double frequency)
        {
            var freqIndex = FindClosestIndexOfFrequency(transform, frequency);
            if (freqIndex < 0 || freqIndex >= transform.FrequencyCount)
                return Enumerable.Empty<double>();
            return transform.EnumerateValuesOfFrequencyAt(freqIndex);
        }

        public static IEnumerable<double> EnumerateValuesOfFrequencyAt(this IBilinearTimeFrequencyRepresentation transform, int index)
        {
            return Enumerable.Range(0, transform.SampleCount).Select(i => transform[i, index]);
        }

        public static IEnumerable<double> EnumerateTimes(this IBilinearTimeFrequencyRepresentation transform)
        {
            return Enumerable.Range(0, transform.SampleCount)
                .Select(i => transform.StartTime + transform.SamplingPeriod * i);
        }

        public static int FindClosestOffsetOfTime(this IBilinearTimeFrequencyRepresentation transform, double time)
        {
            return (int) Math.Round((time - transform.StartTime) / transform.SamplingPeriod);
        }

        public static IEnumerable<double> EnumerateValuesOfTime(this IBilinearTimeFrequencyRepresentation transform, double time)
        {
            var offset = transform.FindClosestOffsetOfTime(time);
            if (offset < 0 || offset >= transform.SampleCount)
                return Enumerable.Empty<double>();
            return transform.EnumerateValuesOfTimeAt(offset);
        }

        public static IEnumerable<double> EnumerateValuesOfTimeAt(this IBilinearTimeFrequencyRepresentation transform, int offset)
        {
            return Enumerable.Range(0, transform.FrequencyCount).Select(i => transform[offset, i]);
        }

        public static IEnumerable<double> EnumerateValues(this IBilinearTimeFrequencyRepresentation transform)
        {
            return from i in Enumerable.Range(0, transform.SampleCount)
                from j in Enumerable.Range(0, transform.FrequencyCount)
                select transform[i, j];
        }
    }
}
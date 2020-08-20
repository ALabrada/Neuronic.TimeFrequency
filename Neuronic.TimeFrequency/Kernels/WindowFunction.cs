using System;

namespace Neuronic.TimeFrequency.Kernels
{
    /// <summary>
    /// Represents a window function.
    /// </summary>
    public class WindowFunction
    {
        private readonly Func<int, int, double> _func;

        internal WindowFunction(Func<int, int, double> func)
        {
            _func = func ?? throw new ArgumentNullException(nameof(func));
        }

        /// <summary>
        /// Evaluates the window function in the specified number of samples.
        /// </summary>
        /// <param name="count">The sample count.</param>
        /// <returns>The window samples.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="count"/> is not positive.</exception>
        public double[] Evaluate(int count)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            var values = new double[count];
            Evaluate(values, 0, count);
            return values;
        }

        /// <summary>
        /// Evaluates the window function in the specified buffer.
        /// </summary>
        /// <param name="values">The buffer.</param>
        /// <param name="start">The start index in <paramref name="values"/>.</param>
        /// <param name="count">The sample count.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> is <c>null</c>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="count"/> is not positive or <paramref name="start"/> is negative.</exception>
        public void Evaluate(double[] values, int start, int count)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            if (start < 0) throw new ArgumentOutOfRangeException(nameof(start));
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));

            for (int i = 0; i < count; i++)
                values[i + start] = _func(i, count - 1);
        }
    }
}
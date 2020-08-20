using System;

namespace Neuronic.TimeFrequency.Kernels
{
    public class WindowFunction
    {
        private readonly Func<int, int, double> _func;

        public WindowFunction(Func<int, int, double> func)
        {
            _func = func ?? throw new ArgumentNullException(nameof(func));
        }

        public double[] Evaluate(int count)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            var values = new double[count];
            Evaluate(values, 0, count);
            return values;
        }

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
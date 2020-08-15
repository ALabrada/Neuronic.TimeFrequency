using System.Collections.Generic;

namespace Neuronic.TimeFrequency
{
    public interface IReadOnlySignal<T> : IReadOnlyList<T>
    {
        double Delay { get; }
        double SamplingRate { get; }
        double SamplingPeriod { get; }
        double Duration { get; }
    }
}
#include "fpproc.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct ProbePoint {
    double x;
    double y;
};

struct LineProbe {
    double x0;
    double y0;
    double x1;
    double y1;
};

struct WindowProbe {
    double xmin;
    double ymin;
    double xmax;
    double ymax;
    int nx;
    int ny;
};

struct CaseSpec {
    std::string circuitName;
    bool hasCircuit = false;
    bool hasPeakWindow = false;
    LineProbe lineProbe{};
    WindowProbe peakWindow{};
    WindowProbe fieldWindow{};
    std::vector<ProbePoint> probes;
};

double bmag(const CMPointVals &u)
{
    return std::sqrt(absq(u.B1) + absq(u.B2));
}

bool getCaseSpec(const std::string &caseName, CaseSpec *spec)
{
    if (caseName == "case1")
    {
        spec->hasCircuit = true;
        spec->circuitName = "coil";
        spec->lineProbe = {61.0, 50.0, 79.0, 50.0};
        spec->fieldWindow = {-5.0, -5.0, 85.0, 105.0, 120, 120};
        spec->probes = {{70.0, 50.0}};
        return true;
    }
    if (caseName == "case2")
    {
        spec->hasCircuit = true;
        spec->circuitName = "coil";
        spec->lineProbe = {40.0, 70.25, 60.0, 70.25};
        spec->hasPeakWindow = true;
        spec->peakWindow = {38.0, 70.2, 62.0, 71.0, 49, 33};
        spec->fieldWindow = {10.0, 0.0, 90.0, 95.0, 120, 120};
        spec->probes = {{50.0, 70.75}, {50.0, 71.5}};
        return true;
    }
    if (caseName == "case3")
    {
        spec->lineProbe = {40.0, 38.0, 40.0, 50.0};
        spec->fieldWindow = {15.0, -5.0, 65.0, 55.0, 120, 120};
        spec->probes = {{40.0, 44.0}, {52.0, 25.0}};
        return true;
    }
    return false;
}

int findCircuitIndex(const FPProc &proc, const std::string &name)
{
    for (int i = 0; i < static_cast<int>(proc.circproplist.size()); ++i)
    {
        if (proc.circproplist[i].CircName == name)
            return i;
    }
    return -1;
}

void writeLineMetrics(FPProc &proc, const LineProbe &lineProbe, std::ofstream &metrics)
{
    proc.contour.clear();
    proc.contour.push_back(CComplex(lineProbe.x0, lineProbe.y0));
    proc.contour.push_back(CComplex(lineProbe.x1, lineProbe.y1));
    CComplex z[4] = {0., 0., 0., 0.};
    proc.LineIntegral(0, z);
    metrics << "line_bn_total," << abs(z[0]) << "\n";
    metrics << "line_bn_avg," << abs(z[1]) << "\n";
}

void writePeakWindowMetric(FPProc &proc, const WindowProbe &window, std::ofstream &metrics)
{
    double peak = 0.;
    for (int iy = 0; iy < window.ny; ++iy)
    {
        for (int ix = 0; ix < window.nx; ++ix)
        {
            const double x = window.xmin + (window.xmax - window.xmin) * ix / (window.nx - 1);
            const double y = window.ymin + (window.ymax - window.ymin) * iy / (window.ny - 1);
            CMPointVals u;
            if (proc.GetPointValues(x, y, u))
                peak = (std::max)(peak, bmag(u));
        }
    }
    metrics << "peak_B_window," << peak << "\n";
}

void writeFieldCsv(FPProc &proc, const WindowProbe &window, std::ofstream &field)
{
    field << "x,y,bmag\n";
    for (int iy = 0; iy < window.ny; ++iy)
    {
        for (int ix = 0; ix < window.nx; ++ix)
        {
            const double x = window.xmin + (window.xmax - window.xmin) * ix / (window.nx - 1);
            const double y = window.ymin + (window.ymax - window.ymin) * iy / (window.ny - 1);
            CMPointVals u;
            double value = 0.;
            if (proc.GetPointValues(x, y, u))
                value = bmag(u);
            field << x << "," << y << "," << value << "\n";
        }
    }
}

} // namespace

int main(int argc, char **argv)
{
    if (argc != 5)
    {
        std::cerr << "Usage: fpproc-vgap <case-name> <ans-path> <metrics-csv> <field-csv>\n";
        return 2;
    }

    const std::string caseName(argv[1]);
    const std::string ansPath(argv[2]);
    const std::string metricsPath(argv[3]);
    const std::string fieldPath(argv[4]);

    CaseSpec spec;
    if (!getCaseSpec(caseName, &spec))
    {
        std::cerr << "Unknown case name: " << caseName << "\n";
        return 2;
    }

    FPProc proc;
    if (!proc.OpenDocument(ansPath))
    {
        std::cerr << "Unable to open ans file: " << ansPath << "\n";
        return 1;
    }

    std::ofstream metrics(metricsPath);
    std::ofstream field(fieldPath);
    if (!metrics || !field)
    {
        std::cerr << "Unable to open metrics or field output file\n";
        return 1;
    }

    metrics << "metric,value\n";

    if (spec.hasCircuit)
    {
        const int circuitIndex = findCircuitIndex(proc, spec.circuitName);
        if (circuitIndex < 0)
        {
            std::cerr << "Circuit not found: " << spec.circuitName << "\n";
            return 1;
        }
        const double amps = abs(proc.circproplist[circuitIndex].Amps);
        const double flux = abs(proc.GetFluxLinkage(circuitIndex));
        metrics << "current," << amps << "\n";
        metrics << "flux_linkage," << flux << "\n";
        if (amps > 0.)
            metrics << "inductance," << (flux / amps) << "\n";
    }

    writeLineMetrics(proc, spec.lineProbe, metrics);

    for (size_t i = 0; i < spec.probes.size(); ++i)
    {
        CMPointVals u;
        double value = 0.;
        if (proc.GetPointValues(spec.probes[i].x, spec.probes[i].y, u))
            value = bmag(u);
        metrics << "probe_" << i << "_Bmag," << value << "\n";
    }

    if (spec.hasPeakWindow)
        writePeakWindowMetric(proc, spec.peakWindow, metrics);

    writeFieldCsv(proc, spec.fieldWindow, field);
    return 0;
}

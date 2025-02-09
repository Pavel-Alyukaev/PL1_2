#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <map>
#include <algorithm>
#include <format>
#include <limits>

struct Point
{
    double x, y, z;

    [[nodiscard]] std::string toString() const
    {
        return std::format("{} {} {}", x, y, z);
    }
};

class Vector
{
public:
    Vector(double x, double y, double z) : m_x(x), m_y(y), m_z(z) {}

    Vector(const Point& a, const Point& b)
            : m_x(b.x - a.x), m_y(b.y - a.y), m_z(b.z - a.z) {}

    Vector operator+(const Vector& other) const {
        return {m_x + other.m_x, m_y + other.m_y, m_z + other.m_z};
    }

    Vector operator*(double scalar) const {
        return {m_x * scalar, m_y * scalar, m_z * scalar};
    }

    [[nodiscard]] double dot(const Vector& other) const {
        return m_x * other.m_x + m_y * other.m_y + m_z * other.m_z;
    }

    [[nodiscard]] Vector cross(const Vector& other) const {
        return {
                m_y * other.m_z - m_z * other.m_y,
                m_z * other.m_x - m_x * other.m_z,
                m_x * other.m_y - m_y * other.m_x
        };
    }

    [[nodiscard]] double length() const {
        return std::sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
    }

    [[nodiscard]] Vector normalized() const {
        double len = length();
        if (len != 0.0) {
            return *this * (1.0 / len);
        }
        return *this;
    }

    [[nodiscard]] Point toPoint() const {
        return {m_x, m_y, m_z};
    }

private:
    double m_x, m_y, m_z;
};

std::vector<Point> readPointsFromFile(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<Point> points;

    if (!file.is_open()) {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << std::endl;
        return points;
    }

    double x, y, z;
    while (file >> x >> y >> z) {
        points.emplace_back(x, y, z);
    }

    file.close();
    return points;
}

enum class stage : int
{
    CheckBluntAngledTriangle,
    MinForBluntAngledTriangle,

};

class TmpData
{
public:
    TmpData() = default;
    ~TmpData() = default;

    void setPoint(const Point& pnt)
    {
        m_point = pnt;
    }

    void addInterval(int i, const Point& pnt1, const Point& pnt2)
    {
        m_intervals.emplace(i, std::make_pair(pnt1, pnt2));
    }

private:
    std::map<int, std::pair<Point, Point>> m_intervals;
    Point m_point;
};





int main() {
    const std::string filename = "C:\\qqqqq\\input.txt";
    std::vector<Point> points = readPointsFromFile(filename);

    if (points.empty()) {
        std::cerr << "Файл пуст или не содержит корректных данных." << std::endl;
        return 1;
    }

    const Point O{2, 0.5, 0.5 };
    std::map<int, Point> closestPoints;
    double minDistance = std::numeric_limits<double>::max();
    const double epsilon = std::numeric_limits<double>::epsilon();

    for (size_t i = 1; i < points.size(); ++i) {
        Vector v1(points[i - 1], O );
        Vector v2(O, points[i]);
        Vector v3(points[i - 1], points[i]);
        auto q = v3.dot(v1);


        double currentMin;
        if(q <= 0.0)
        {
            currentMin = v1.length();
        }
        else if( q >= v3.length())
        {
            currentMin = v1.length();
        }
        else
            closestPoints.emplace(i, (v3 * ( q / v3.length()) + v1).toPoint());

        if (currentMin - minDistance > epsilon) {
            continue;
        } else if (minDistance - currentMin > epsilon) {
            minDistance = currentMin;
            closestPoints.clear();
        }






//
//
//
//
//
//        double lengths[] = {v1.length(), v2.length(), v3.length() / 2.0};
//        auto minIt = std::min_element(std::begin(lengths), std::end(lengths));
//        double currentMin = *minIt;
//
//        if (currentMin - minDistance > epsilon) {
//            continue;
//        } else if (minDistance - currentMin > epsilon) {
//            minDistance = currentMin;
//            closestPoints.clear();
//        }
//
//        size_t index = std::distance(std::begin(lengths), minIt);
//        if (index == 0) {
//            closestPoints.emplace(i, points[i - 1]);
//        } else if (index == 1) {
//            closestPoints.emplace(i, points[i]);
//        } else {
//            closestPoints.emplace(i, (v3 * 0.5 + Vector({0, 0, 0}, O)).toPoint());
//        }
    }

    for (const auto& [index, point] : closestPoints) {
        std::cout << "segment " << index << " point " << point.toString() << std::endl;
    }

    return 0;
}
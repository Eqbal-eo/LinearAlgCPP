#include <iostream>
#include <vector>

class Imaginary
{
private:
    double imagPart;

public:
    Imaginary(double imag = 0.0) : imagPart(imag) {}

    double getImagPart() const
    {
        return imagPart;
    }

    void setImagPart(double imag)
    {
        imagPart = imag;
    }

    void display() const
    {
        std::cout << imagPart << "i";
    }
};

class Complex : public Imaginary
{
private:
    double real;

public:
    Complex(double real = 0.0, double imag = 0.0) : real(real), Imaginary(imag) {}

    Complex operator+(const Complex &other) const
    {
        return Complex(real + other.real, Imaginary::getImagPart() + other.getImagPart());
    }

    Complex operator-(const Complex &other) const
    {
        return Complex(real - other.real, Imaginary::getImagPart() - other.getImagPart());
    }

    Complex operator*(const Complex &other) const
    {
        return Complex(real * other.real - Imaginary::getImagPart() * other.getImagPart(),
                       real * other.getImagPart() + Imaginary::getImagPart() * other.real);
    }

    Complex operator/(const Complex &other) const
    {
        double denominator = other.real * other.real + other.getImagPart() * other.getImagPart();
        return Complex((real * other.real + Imaginary::getImagPart() * other.getImagPart()) / denominator,
                       (Imaginary::getImagPart() * other.real - real * other.getImagPart()) / denominator);
    }

    Complex operator+=(const Complex &other)
    {
        real += other.real;
        Imaginary::setImagPart(Imaginary::getImagPart() + other.getImagPart());
        return *this;
    }

    Complex operator-=(const Complex &other)
    {
        real -= other.real;
        Imaginary::setImagPart(Imaginary::getImagPart() - other.getImagPart());
        return *this;
    }

    Complex operator*=(const Complex &other)
    {
        Complex temp = *this * other;
        real = temp.real;
        Imaginary::setImagPart(temp.getImagPart());
        return *this;
    }

    Complex operator/=(const Complex &other)
    {
        Complex temp = *this / other;
        real = temp.real;
        Imaginary::setImagPart(temp.getImagPart());
        return *this;
    }

    bool operator==(const Complex &other) const
    {
        return (real == other.real) && (Imaginary::getImagPart() == other.getImagPart());
    }

    void display() const
    {
        std::cout << real << " + ";
        Imaginary::display();
        std::cout << " ";
    }
};

class Vector
{
private:
    size_t size;
    std::vector<Complex> elements;

public:
    Vector(size_t size = 0, Complex element = Complex(0.0, 0.0)) : size(size), elements(size, element) {}

    Vector operator+(const Vector &other) const
    {
        Vector result(size);
        for (size_t i = 0; i < size; ++i)
        {
            result.elements[i] = elements[i] + other.elements[i];
        }
        return result;
    }

    Vector operator-(const Vector &other) const
    {
        Vector result(size);
        for (size_t i = 0; i < size; ++i)
        {
            result.elements[i] = elements[i] - other.elements[i];
        }
        return result;
    }

    Complex operator*(const Vector &other) const
    {
        Complex result(0.0, 0.0);
        for (size_t i = 0; i < size; ++i)
        {
            result += elements[i] * other.elements[i];
        }
        return result;
    }

    Vector &operator=(const Vector &other)
    {
        if (this != &other)
        {
            size = other.size;
            elements = other.elements;
        }
        return *this;
    }

    Vector &operator+=(const Vector &other)
    {
        for (size_t i = 0; i < size; ++i)
        {
            elements[i] += other.elements[i];
        }
        return *this;
    }

    Vector &operator-=(const Vector &other)
    {
        for (size_t i = 0; i < size; ++i)
        {
            elements[i] -= other.elements[i];
        }
        return *this;
    }

    bool operator==(const Vector &other) const
    {
        return size == other.size && elements == other.elements;
    }

    void display() const
    {
        std::cout << "[ ";
        for (const auto &element : elements)
        {
            element.display();
        }
        std::cout << "]" << std::endl;
    }
};

class Matrix
{
protected:
    size_t row;
    size_t col;
    std::vector<std::vector<Complex>> elements;

public:
    Matrix(size_t row = 0, size_t col = 0, Complex element = Complex(0.0, 0.0))
        : row(row), col(col), elements(row, std::vector<Complex>(col, element)) {}

    Matrix operator+(const Matrix &other) const
    {
        Matrix result(row, col);
        for (size_t i = 0; i < row; ++i)
        {
            for (size_t j = 0; j < col; ++j)
            {
                result.elements[i][j] = elements[i][j] + other.elements[i][j];
            }
        }
        return result;
    }

    Matrix operator-(const Matrix &other) const
    {
        Matrix result(row, col);
        for (size_t i = 0; i < row; ++i)
        {
            for (size_t j = 0; j < col; ++j)
            {
                result.elements[i][j] = elements[i][j] - other.elements[i][j];
            }
        }
        return result;
    }

    Matrix operator*(const Matrix &other) const
    {
        if (col != other.row)
        {
            std::cerr << "Matris boyutları uyumsuz!" << std::endl;
            return Matrix();
        }

        Matrix result(row, other.col);
        for (size_t i = 0; i < row; ++i)
        {
            for (size_t j = 0; j < other.col; ++j)
            {
                for (size_t k = 0; k < col; ++k)
                {
                    result.elements[i][j] += elements[i][k] * other.elements[k][j];
                }
            }
        }
        return result;
    }

    Matrix &operator=(const Matrix &other)
    {
        if (this != &other)
        {
            row = other.row;
            col = other.col;
            elements = other.elements;
        }
        return *this;
    }

    Matrix &operator+=(const Matrix &other)
    {
        for (size_t i = 0; i < row; ++i)
        {
            for (size_t j = 0; j < col; ++j)
            {
                elements[i][j] += other.elements[i][j];
            }
        }
        return *this;
    }

    Matrix &operator-=(const Matrix &other)
    {
        for (size_t i = 0; i < row; ++i)
        {
            for (size_t j = 0; j < col; ++j)
            {
                elements[i][j] -= other.elements[i][j];
            }
        }
        return *this;
    }

    bool operator==(const Matrix &other) const
    {
        return row == other.row && col == other.col && elements == other.elements;
    }

    void display() const
    {
        for (size_t i = 0; i < row; ++i)
        {
            std::cout << "[ ";
            for (size_t j = 0; j < col; ++j)
            {
                elements[i][j].display();
            }
            std::cout << "]" << std::endl;
        }
        std::cout << std::endl;
    }

    Matrix transpose() const
    {
        Matrix result(col, row);
        for (size_t i = 0; i < row; ++i)
        {
            for (size_t j = 0; j < col; ++j)
            {
                result.elements[j][i] = elements[i][j];
            }
        }
        return result;
    }

    Complex determinant() const
    {
        if (row != col)
        {
            std::cerr << "Determinant hesaplamak için kare matris gerekli!" << std::endl;
            return Complex();
        }

        if (row == 1)
        {
            return elements[0][0];
        }
        else if (row == 2)
        {
            return elements[0][0] * elements[1][1] - elements[0][1] * elements[1][0];
        }
        else
        {
            Complex det(0.0, 0.0);
            for (size_t i = 0; i < col; ++i)
            {
                Matrix submatrix(row - 1, col - 1);
                for (size_t j = 1; j < row; ++j)
                {
                    for (size_t k = 0, l = 0; k < col; ++k)
                    {
                        if (k == i)
                            continue;
                        submatrix.elements[j - 1][l++] = elements[j][k];
                    }
                }
                Complex sign = (i % 2 == 0) ? Complex(1.0, 0.0) : Complex(-1.0, 0.0);
                det = det + (sign * elements[0][i] * submatrix.determinant());
            }
            return det;
        }
    }
};

class SquareMatrix : public Matrix
{
public:
    SquareMatrix(size_t size = 0, Complex element = Complex(0.0, 0.0)) : Matrix(size, size, element) {}
};

class IdentityMatrix : public SquareMatrix
{
public:
    IdentityMatrix(size_t size = 0) : SquareMatrix(size)
    {
        for (size_t i = 0; i < size; ++i)
        {
            elements[i][i] = Complex(1.0, 0.0);
        }
    }
};

class TriangleMatrix : public SquareMatrix
{
public:
    TriangleMatrix(size_t size = 0, Complex element = Complex(0.0, 0.0)) : SquareMatrix(size, element)
    {
        for (size_t i = 0; i < size; ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                elements[i][j] = Complex(0.0, 0.0);
            }
        }
    }
};

class LinearAlgebraObject
{
public:
    virtual void print() const = 0;
    virtual ~LinearAlgebraObject() = default;
};

Imaginary::Imaginary(double imag) : imagPart(imag) {}

Complex::Complex(double real, double imag) : real(real), Imaginary(imag) {}

Vector::Vector(size_t size, Complex element) : size(size), elements(size, element) {}

Vector Vector::operator+(const Vector &other) const
{
    Vector result(size);
    for (size_t i = 0; i < size; ++i)
    {
        result.elements[i] = elements[i] + other.elements[i];
    }
    return result;
}

Vector Vector::operator-(const Vector &other) const
{
    Vector result(size);
    for (size_t i = 0; i < size; ++i)
    {
        result.elements[i] = elements[i] - other.elements[i];
    }
    return result;
}

Complex Vector::operator*(const Vector &other) const
{
    Complex result(0.0, 0.0);
    for (size_t i = 0; i < size; ++i)
    {
        result += elements[i] * other.elements[i];
    }
    return result;
}

Vector &Vector::operator=(const Vector &other)
{
    if (this != &other)
    {
        size = other.size;
        elements = other.elements;
    }
    return *this;
}

Vector &Vector::operator+=(const Vector &other)
{
    for (size_t i = 0; i < size; ++i)
    {
        elements[i] += other.elements[i];
    }
    return *this;
}

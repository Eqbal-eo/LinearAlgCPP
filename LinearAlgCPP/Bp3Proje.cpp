
/**
* @file Bp3 Proje
* @lineer cebir kütüphanesi oluşturmayı amaçlamaktadır. 
* @2. Odev
* @date 29/12/2023
* @author Mhd Eqbal AYOUB/ mhdeqbal.ayoub@stu.fsm.edu.tr
*/

#include <iostream>
#include <vector>

using namespace std;

// LinearAlgebraObject soyut sınıfı
class LinearAlgebraObject
{
public:
    // Kurucu (constructor)
    LinearAlgebraObject() {}

    LinearAlgebraObject(const LinearAlgebraObject &other) {}

    virtual void print() const = 0;

    virtual ~LinearAlgebraObject() {}
};

// Imaginary sınıfı
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

    void print() const
    {
        cout << imagPart << "i";
    }

    Imaginary operator+(const Imaginary &other) const
    {
        return Imaginary(imagPart + other.imagPart);
    }

    Imaginary operator-(const Imaginary &other) const
    {
        return Imaginary(imagPart - other.imagPart);
    }

    Imaginary operator*(const Imaginary &other) const
    {
        return Imaginary(imagPart * other.imagPart);
    }

    bool operator==(const Imaginary &other) const
    {
        return imagPart == other.imagPart;
    }
};

// Complex sınıfı
class Complex : Imaginary
{
private:
    double real;
    double imaginary;

public:
    Complex(double real = 0.0, double imaginary = 0.0)
        : real(real), imaginary(imaginary) {}

    Complex operator+(const Complex &other) const
    {
        return Complex(real + other.real, imaginary + other.imaginary);
    }

    Complex operator-(const Complex &other) const
    {
        return Complex(real - other.real, imaginary - other.imaginary);
    }

    Complex operator*(const Complex &other) const
    {
        return Complex(real * other.real - imaginary * other.imaginary,
                       real * other.imaginary + imaginary * other.real);
    }

    bool operator==(const Complex &other) const
    {
        return (real == other.real) && (imaginary == other.imaginary);
    }

    void print() const
    {
        std::cout << real << " + " << imaginary << "i"
                  << " ";
    }
};


// Vectr sınıfı
class Vector
{
private:
    size_t size;
    vector<Complex> elements;

public:
    Vector(size_t size = 0, Complex element = Complex(0.0, 0.0))
        : size(size), elements(size, element) {}

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
            result = result + (elements[i] * other.elements[i]);
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
            elements[i] = elements[i] + other.elements[i];
        }
        return *this;
    }

    Vector &operator-=(const Vector &other)
    {
        for (size_t i = 0; i < size; ++i)
        {
            elements[i] = elements[i] - other.elements[i];
        }
        return *this;
    }

    bool operator==(const Vector &other) const
    {
        return size == other.size && elements == other.elements;
    }

    void print() const
    {
        cout << "[ ";
        for (const auto &element : elements)
        {
            element.print();
        }
        cout << "]" << endl;
    }
};

// Matric sınıfı
class Matrix
{
protected:
    size_t row;
    size_t col;
    vector<vector<Complex>> elements;

public:
    Matrix(size_t row = 0, size_t col = 0, Complex element = Complex(0.0, 0.0))
        : row(row), col(col), elements(row, vector<Complex>(col, element)) {}

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
            cerr << "Matris boyutları uyumsuz!" << endl;
            return Matrix();
        }

        Matrix result(row, other.col);
        for (size_t i = 0; i < row; ++i)
        {
            for (size_t j = 0; j < other.col; ++j)
            {
                for (size_t k = 0; k < col; ++k)
                {
                    result.elements[i][j] = result.elements[i][j] + (elements[i][k] * other.elements[k][j]);
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
                elements[i][j] = elements[i][j] + other.elements[i][j];
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
                elements[i][j] = elements[i][j] - other.elements[i][j];
            }
        }
        return *this;
    }

    bool operator==(const Matrix &other) const
    {
        return row == other.row && col == other.col && elements == other.elements;
    }

    void print() const
    {
        for (size_t i = 0; i < row; ++i)
        {
            cout << "[ ";
            for (size_t j = 0; j < col; ++j)
            {
                elements[i][j].print();
            }
            cout << "]" << endl;
        }
        cout << endl;
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
            cerr << "Determinant hesaplamak için kare matris gerekli!" << endl;
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

// SquareMatrix sınıfı
class SquareMatrix : public Matrix
{
public:
    SquareMatrix(size_t size = 0, Complex element = Complex(0.0, 0.0))
        : Matrix(size, size, element) {}
};

// IdentityMatrix sınıfı
class IdentityMatrix : public SquareMatrix
{
public:
    IdentityMatrix(size_t size = 0)
        : SquareMatrix(size)
    {
        for (size_t i = 0; i < size; ++i)
        {
            elements[i][i] = Complex(1.0, 0.0);
        }
    }
};

// TriangleMatrix sınıfı
class TriangleMatrix : public SquareMatrix
{
public:
    TriangleMatrix(size_t size = 0, Complex element = Complex(0.0, 0.0))
        : SquareMatrix(size, element)
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

int main()
{
    // Vector örneği
    Vector vector1(3, Complex(1.0, 2.0));
    Vector vector2(3, Complex(3.0, 4.0));

    Vector sumVector = vector1 + vector2;
    cout << "\n" <<"Vector Toplama: "
         << "\n";
    sumVector.print();

    Vector diffVector = vector1 - vector2;
    cout << "\n"
         << "Vector Cikarma: "
         << "\n";
    diffVector.print();

    Complex dotProduct = vector1 * vector2;
    cout << "\n"
         << "Vector Carpim: "
         << "\n";
    dotProduct.print();

    Vector assignedVector = vector1;
    assignedVector += vector2;
    cout << "\n"
         << "Vector Atama ve Toplama Atama: "
         << "\n";
    assignedVector.print();

    // Matrix örneği
    Matrix matrix1(2, 2, Complex(1.0, 2.0));
    Matrix matrix2(2, 2, Complex(3.0, 4.0));

    Matrix sumMatrix = matrix1 + matrix2;
    cout << "\n"
         << "Matrix Toplama: "
         << "\n";
    sumMatrix.print();

    Matrix diffMatrix = matrix1 - matrix2;
    cout << "\n"
         << "Matrix Cikarma: "
         << "\n";
    diffMatrix.print();

    Matrix productMatrix = matrix1 * matrix2;
    cout << "\n"
         << "Matrix Carpma: "
         << "\n";
    productMatrix.print();

    Matrix assignedMatrix = matrix1;
    assignedMatrix += matrix2;
    cout << "\n"
         << "Matrix Atama ve Toplama Atama: "
         << "\n";
    assignedMatrix.print();

    Matrix transposedMatrix = matrix1.transpose();
    cout << "\n"
         << "Matrix Transpozu: "
         << "\n";
    transposedMatrix.print();

    Complex determinantValue = matrix1.determinant();
    cout << "\n"
         << "Matrix Determinanti: "
         << "\n";
    determinantValue.print();

    // SquareMatrix örneği
    SquareMatrix squareMatrix(3, Complex(1.0, 2.0));
    cout << "SquareMatrix:" << endl;
    squareMatrix.print();

    // IdentityMatrix örneği
    IdentityMatrix identityMatrix(3);
    cout << "IdentityMatrix:" << endl;
    identityMatrix.print();

    // TriangleMatrix örneği
    TriangleMatrix triangleMatrix(3, Complex(1.0, 2.0));
    cout << "TriangleMatrix:" << endl;
    triangleMatrix.print();

    return 0;
}
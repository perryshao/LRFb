#ifndef Unit2H
#define Unit2H

typedef long Z;

Z max(Z a, Z b)
{
    if (a > b)
        return a;
    else
        return b;
}
Z min(Z a, Z b)
{
    if (a < b)
        return a;
    else
        return b;
}

//int max(int a, int b) { if ( a > b  ) return a; else  return  b; }
//int min(int a, int b) { if ( a < b  ) return a; else  return  b; }

class Q
{
  public:
    Z n;
    Z d;
    Q(int a, int b)
    {
        n = a;
        d = b;
    }
    Q(int a)
    {
        n = a;
        d = 1;
    }
    Q() {};
    Q(const Q &q)
    {
        n = q.n;
        d = q.d;
    }
    Q &operator=(const Q &q)
    {
        n = q.n;
        d = q.d;
        return *this;
    }
};

bool operator>(const Q &a, const Q &b)
{
    return (a.n * b.d > b.n * a.d);
}
bool operator<(const Q &a, const Q &b)
{
    return (a.n * b.d < b.n * a.d);
}
bool operator>=(const Q &a, const Q &b)
{
    return (a.n * b.d >= b.n * a.d);
}
bool operator<=(const Q &a, const Q &b)
{
    return (a.n * b.d <= b.n * a.d);
}
bool operator==(const Q &a, const Q &b)
{
    return (a.n * b.d == b.n * a.d);
}
bool operator!=(const Q &a, const Q &b)
{
    return (a.n * b.d != b.n * a.d);
}

class point
{
  public:
    Z x;
    Z y;

    point(int X, int Y) : x(X), y(Y) {}
    point() {}
    point(const point &p) : x(p.x), y(p.y) {}
    point &operator=(const point &p)
    {
        x = p.x;
        y = p.y;
        return *this;
    }
    void prime()
    {
        int n = min(abs(x), abs(y)); // reduction gcd=1  // to replace
        for (int d = 2; d <= n; d++)
            while ((x % d == 0) && (y % d == 0))
            {
                x /= d;
                y /= d;
            }
    }
};

typedef point vector;

point operator-(const point &a, const point &b)
{
    point r;
    r.x = a.x - b.x;
    r.y = a.y - b.y;
    return r;
}

point operator+(const point &a, const point &b)
{
    point r;
    r.x = a.x + b.x;
    r.y = a.y + b.y;
    return r;
}

// vector product

int operator^(const point &a, const point &b)
{
    return a.x * b.y - a.y * b.x;
}

// dot product

int operator*(const point &a, const point &b)
{
    return a.x * b.x + a.y * b.y;
}

typedef point *Lpoint;

//---------------------------------------------------------------------------
#endif

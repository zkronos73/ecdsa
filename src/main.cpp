#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <gmp.h>
#include <gmpxx.h>
#include <memory>
#include <stdexcept>
#include <nlohmann/json.hpp>

#include <alt_bn128.hpp>

using namespace AltBn128;
using json = nlohmann::json;

#define handle_error(msg) \
           do { perror(msg); exit(EXIT_FAILURE); } while (0)


class mpz_tools
{
    public:
        static mpz_class inverse(const mpz_class &op1, const mpz_class &op2) { 
            mpz_class dst;
            mpz_invert(dst.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t()); 
            return dst;
            };
        static G1Point mul(G1Point &p, const mpz_class &scalar) {
            uint8_t _scalar[32];
            G1Point result;
            for (int i=0;i<32;i++) _scalar[i] = 0;
            mpz_export((void *)_scalar, NULL, -1, 8, -1, 0, scalar.get_mpz_t());
            G1.mulByScalar(result, p, _scalar, 32);
            return result;
        }
        static G1Point makePoint(const mpz_class &x, const mpz_class &y) {
            G1PointAffine a;
            G1Point result;
            G1.F.fromMpz(a.x, x.get_mpz_t());
            G1.F.fromMpz(a.y, y.get_mpz_t());
            G1.copy(result, a);
            return result;
        }
        static mpz_class get(const F1Element &element) {
            mpz_class result;
            G1.F.toMpz(result.get_mpz_t(), element);
            return result;
        }
        static F1Element toField(const mpz_class &a) {
            F1Element result;
            G1.F.fromMpz(result, a.get_mpz_t());
            return result;
        }
        static mpz_class getX(G1Point &point) {
            G1PointAffine pa;
            G1.copy(pa, point);
            return get(pa.x);
        }
        static mpz_class getY(G1Point &point) {
            G1PointAffine pa;
            G1.copy(pa, point);
            return get(pa.y);
        }
        static mpz_class sqrt(const mpz_class &op1) {
            mpz_class result;
            mpz_sqrt(result.get_mpz_t(), op1.get_mpz_t());
            return result;
        }
};

/*

long pow_mod(long x, long n, long p) {
  if (n == 0) return 1;
  if (n & 1)
    return (pow_mod(x, n-1, p) * x) % p;
  x = pow_mod(x, n/2, p);
  return (x * x) % p;
}
*/

/* Takes as input an odd prime p and n < p and returns r
 * such that r * r = n [mod p]. */
/*void tonelliShanks ( F1Element &r, const F1Element &n, const F1Element &p )
 {
    F1Element s = G1.F.zero();
    F1Element q = G1.F.sub(G1.F.zero(), G1.F.one());

    while ((q & 1)) {

    }

}*/

/**
 * https://gist.github.com/LaurentMazare/6745649
 *
 * */

mpz_class pow ( const mpz_class &x, const mpz_class &n, const mpz_class &p ) 
{
    if (n == 0) return 1;
    if ((n & 1) == 1) {
        return (pow(x, n-1, p) * x) % p;
    }
    mpz_class x2 = pow(x, n/2, p);
    return (x2 * x2) % p;
}

mpz_class sqrtTonelliShanks ( const mpz_class &n, const mpz_class &p )
{
    mpz_class s = 1;
    mpz_class q = p - 1;
    while ((q & 1) == 0) {
        q = q / 2;
        ++s;
    }
    if (s == 1) {
        mpz_class r = pow(n, (p+1)/4, p);
        if ((r * r) % p == n) return r;
        return 0;
    }

    mpz_class z = 1;
    while (pow(++z, (p - 1)/2, p) != (p - 1));
//    std::cout << "Z found: " << z << "\n";
    mpz_class c = pow(z, q, p);
    mpz_class r = pow(n, (q+1)/2, p);
    mpz_class t = pow(n, q, p);
    mpz_class m = s;
    while (t != 1) {
        mpz_class tt = t;
        mpz_class i = 0;
        while (tt != 1) {
            tt = (tt * tt) % p;
            ++i;
            if (i == m) return 0;
        }
        mpz_class b = pow(c, pow(2, m-i-1, p-1), p);
        mpz_class b2 = (b * b) % p;
        r = (r * b) % p;
        t = (t * b2) % p;
        c = b2;
        m = i;
    }
    if (((r * r) % p) == n) return r;
    return 0;
}

int main(int argc, char **argv) 
{
    mpz_class n("21888242871839275222246405745257275088548364400416034343698204186575808495617");
    mpz_class q("21888242871839275222246405745257275088696311157297823662689037894645226208583");

    mpz_class secret("21888242871839275222246405745250000088548364400416034343698204186575808495617");
    mpz_class hash("234566547655868674654654");

    /** sign **/
    mpz_class r(0), s(0);
    G1Point randomPoint;
    G1Point publicPoint = mpz_tools::mul(G1.one(), secret);
    std::cout << "PUBLIC: " << G1.toString(publicPoint) << "\n";

    gmp_randclass rnd(gmp_randinit_default);
    while (r == 0 || s == 0) {
        mpz_class rnum = rnd.get_z_range(n);
        // rnum = 1;
        std::cout << "RND:" << rnum << "\n";

        randomPoint = mpz_tools::mul(G1.one(), rnum);
        std::cout << "randomPoint: " << G1.toString(randomPoint) << "\n";
        std::cout << "X:" << mpz_tools::getX(randomPoint) << "\n";
        r = mpz_tools::getX(randomPoint) % n;        
        // G1.F.toMpz(r.get_mpz_t(), randomPoint.x);
        // r = r % n;
        std::cout << "r=" << r << "\n";

        s = (((hash + r * secret)) * mpz_tools::inverse(rnum, n)) % n;
        std::cout << "s=" << s << "\n";
    }
    std::cout << "X:" << mpz_tools::getX(randomPoint) << " Y:" << mpz_tools::getY(randomPoint) << "\n";
    mpz_class y = mpz_tools::getY(randomPoint);
    // G1.F.toMpz(y.get_mpz_t(), randomPoint.y);
    // hash = "234566547555868674654654";

    std::cout << "\n";
    mpz_class recoverId = y & 1;
    if (mpz_tools::getX(randomPoint) >= n) {
        recoverId += 2;
    }
    std::cout << "r:" << r << " s:" << s << " recoverId:" << recoverId << "\n";


    bool result = true;
    while (result) {

        if (r < 1 || r > (n - 1)) {
            result = false;
            std::cout << "Fail on R verification\n";
            break;
        }

        if (s < 1 || s > (n - 1)) {
            result = false;
            std::cout << "Fail on s verification\n";
            break;
        }

        mpz_class inv_r = mpz_tools::inverse(r, n);
        mpz_class x = r;
        // x = 1;
        if ((recoverId & 2) != 0) {
            x += n;
        }

        // y² = x³ + 7 (secp256k1)    
        // y² = x³ + 3 (bn128)

        std::cout << "[x: " << x << "]\n";

        mpz_class yy = (((x * x)* x) + 3) % q;
        mpz_class _y = sqrtTonelliShanks(yy, q);
        if ((recoverId & 1) != (_y & 1)) {
            _y = q - _y;
        }

        std::cout << "[y:" << _y << "]\n";
  

        printf("**** METHOD 2 x MUL + SUB *****");

        G1Point x1y1 = mpz_tools::makePoint(x, _y);
        std::cout << "xy:" << G1.toString(x1y1) << "\n";
        G1Point u1 = mpz_tools::mul(G1.one(), (hash * inv_r) % n);
        G1Point u2 = mpz_tools::mul(x1y1, (s * inv_r) % n);
        G1Point v;
        G1.sub(v, u2, u1);
        std::cout << "result:" << G1.toString(v) << "\n";

        // FAKE to verify others scalars
        /* u1 = mpz_tools::mul(G1.one(), (s * inv_r) % n);
        u2 = mpz_tools::mul(x1y1, (hash * inv_r) % n);
        G1.sub(v, u2, u1);
        std::cout << "fake-result:" << G1.toString(v) << "\n";
*/
        printf("**** METHOD  *****");
        G1Point _x1y1;
//        G1Point _x1y1 = mpz_tools::makePoint(x, q - _y);

        std::cout << "xy:" << G1.toString(x1y1) << "\n";

        /*
        G1Point u1 = mpz_tools::mul(G1.one(), (hash * inv_r) % n);
        G1Point u2 = mpz_tools::mul(x1y1, (s * inv_r) % n);
        G1Point v;
        G1.sub(v, u2, u1);
        */
        std::cout << "p0:" << G1.toString(G1.zero()) << "\n";

        G1Point g1_x1y1;
        G1.neg(_x1y1, x1y1);
        G1.add(g1_x1y1, G1.one(), _x1y1);

        std::cout << "p1:" << G1.toString(G1.one()) << "\n";
        std::cout << "p2:" << G1.toString(_x1y1) << "\n";
        std::cout << "p12:" << G1.toString(g1_x1y1) << "\n";
        

//        mpz_class s1 = (hash * inv_r) % n;
//        mpz_class s2 = (s * inv_r) % n;
        mpz_class s1 = (hash * inv_r) % 10000;
        mpz_class s2 = (s * inv_r) % 10000;
        s1 = 2;
        s2 = 3;
        std::cout << "k1: " << s1 << "\nk2: " << s2 << "\n";

        G1Point res = G1.zero();
        for (int bit = 16; bit >= 0; --bit) {
            int b1 = mpz_tstbit(s1.get_mpz_t(), bit);
            int b2 = mpz_tstbit(s2.get_mpz_t(), bit);
            std::cout << "#A bit:" << bit << " b1:" << b1 << " b2:" << b2 << "\n";
            if (b1 && b2) {
                std::cout << "#B 1+1 res:" << G1.toString(res) << "\n";
                std::cout << "#B 1+1 g1_x1y1:" << G1.toString(g1_x1y1) << "\n";
                G1.add(res, res, g1_x1y1);
            }
            else if (b1) {
                std::cout << "#B 1+0 res:" << G1.toString(res) << "\n";
                std::cout << "#B 1+0 g1_x1y1:" << G1.toString(G1.one()) << "\n";
                G1.add(res, res, G1.one());
            }
            else if (b2) {
                std::cout << "#B 0+1 res:" << G1.toString(res) << "\n";
                std::cout << "#B 0+1 g1_x1y1:" << G1.toString(_x1y1) << "\n";
                G1.add(res, res, _x1y1);
            }
            std::cout << "#C res:" << G1.toString(res) << "\n";
            if (bit) G1.dbl(res, res);
            std::cout << "#D res:" << G1.toString(res) << "\n";
        }
        std::cout << "xy(res):" << G1.toString(res) << "\n";
        G1.neg(res, res);

        std::cout << "\n";
/*        G1Point u1 = mpz_tools::mul(G1.one(), (hash * inv_r) % n);
        G1Point u2 = mpz_tools::mul(x1y1, (s * inv_r) % n);
        G1Point v;
        G1.sub(v, u2, u1);*/


        exit(EXIT_FAILURE);
/*        
        if (mpz_tools::getY(v) == 0) {
            result = false;
            std::cout << "y == 0\n";
            break;
        }

        if ((mpz_tools::getX(v) % n) != r) {
            std::cout << "x[" << mpz_tools::getX(v) << "] % n[" << n << "] != r[" << r << "]\n";
            result = false;
            break;
        }
        break;*/
    }
//    std::cout << "RESULT:" << result << "\n";


/*
    bool result = true;
    while (result) {

        if (r < 1 || r > (n - 1)) {
            result = false;
            std::cout << "Fail on R verification\n";
            break;
        }

        if (s < 1 || s > (n - 1)) {
            result = false;
            std::cout << "Fail on s verification\n";
            break;
        }

        mpz_class inv = mpz_tools::inverse(s, n);
        G1Point u1 = mpz_tools::mul(G1.one(), (hash * inv) % n);
        G1Point u2 = mpz_tools::mul(publicPoint, (r * inv) % n);
        G1Point v;
        G1.add(v, u1, u2);
        if (mpz_tools::getY(v) == 0) {
            result = false;
            std::cout << "y == 0\n";
            break;
        }

        if ((mpz_tools::getX(v) % n) != r) {
            std::cout << "x[" << mpz_tools::getX(v) << "] % n[" << n << "] != r[" << r << "]\n";
            result = false;
            break;
        }
        break;
    }
    std::cout << "RESULT:" << result << "\n";
*/
//    G1.mulByScalar()

        // mpz_set_str(altBbn128r, "3", 10);

    // mpz_init(secret);
    // mpz_set_str(secret, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

    // mpz_sub_ui(curveOrder, altBbn128r, 1);
    // uint64_t value[4] = { 0x0FFFFFFFFFFFFFFF, 0x0FFFFFFFFFFFFFFF, 0x0FFFFFFFFFFFFFFF, 0x0FFFFFFFFFFFFFFF };
    
    exit(EXIT_SUCCESS);
}

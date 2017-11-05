//
//  main.c
//  Shad Sterling - MATH 4620 - Project 2
//
//  Created by Shad Sterling on 2014-Apr-23.
//  Copyright (c) 2014 Shad Sterling. All rights reserved.
//

/*
 
 Tested & Working in Xcode 5.1.1
 
 Gives unicode output (may not render properly on old terminals)
 
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

typedef long double ld;
typedef unsigned int ui;

ui udiff( ui i, ui j ) { return i>j ? i-j : j-i; }
ui umax( ui f, ui g ) { return f>g ? f : g; }
ui umin( ui f, ui g ) { return f<g ? f : g; }
ui usub( ui f, ui g ) { return f>g ? f-g : 0; } // subtraction with floor

void crash() { *(int*)1 = *(int*)1; } // this will make it crash, right?

struct matrix_s {
	ld (*get)( struct matrix_s *m, ui i, ui j ); //= &matrix_get_zero;
	bool (*put)( struct matrix_s *m, ui i, ui j, ld a ); //=&matrix_put_fail;
	ui size; //=n
	void *data; //ld**
}; typedef struct matrix_s matrix;
ld matrix_get_zero( matrix *m, ui i, ui j ) { return 0; }
bool matrix_put_fail( matrix *m, ui i, ui j, ld a ) { return false; }

struct vector_s {
	ld (*get)( struct vector_s *v, ui i ); //= &vector_get_zero;
	bool (*put)( struct vector_s *v, ui i, ld a ); //=&vector_put_fail;
	ui size; //=n
	void *data; //ld*
}; typedef struct vector_s vector;
ld vector_get_zero( vector *v, ui i ) { return 0; }
bool vector_put_fail( vector *v, ui i, ld a ) { return false; }

int main( int argc, const char * argv[] );
ld elemA( ui n, ui i, ui j );
ld elemB( ui n, ui i );

vector *solve_crout( matrix *A, vector *b );
matrix **factor_crout( matrix *A, ui k );
vector *fwdsub( matrix *L, vector *b );
vector *baksub( matrix *U, vector *b );

struct gauss_info { vector *x, *b, *r; ld n; ui i; };
struct gauss_info solve_gauss( matrix *A, vector *b, vector *x0, ld tolerance, ui ceiling );
vector *gauss_iterate( vector *xn, matrix *A, vector *b, vector *xp );
struct gauss_info_fast { ld n; ui i; };
struct gauss_info_fast solve_gauss_fast( matrix *A, vector *b, vector *x0, ld tolerance, ui ceiling );

ld matrix_get_elemA( matrix *m, ui i, ui j ) { return elemA( m->size, i, j ); }
void matrix_print( matrix *m );

vector *matrix_apply_vector_col( matrix *A, vector *x );
ld matrix_norm_inf( matrix *A );

matrix *matrix_new_full( ui n );
ld matrix_get_full( matrix *m, ui i, ui j );
bool matrix_put_full( matrix *m, ui i, ui j, ld a );
void matrix_free_full( matrix *m );

struct matrix_data_diag { ui k, s; ld* e; };
matrix *matrix_new_diag( ui n, ui k );
void matrix_free_diag( matrix *m );

matrix *matrix_new_lower( ui n, ui k );
ld matrix_get_lower( matrix *m, ui i, ui j );
bool matrix_put_lower( matrix *m, ui i, ui j, ld a );
void matrix_free_lower( matrix *m ) { matrix_free_diag(m); }

matrix *matrix_new_upper( ui n, ui k );
ld matrix_get_upper( matrix *m, ui i, ui j );
bool matrix_put_upper( matrix *m, ui i, ui j, ld a );
void matrix_free_upper( matrix *m ) { matrix_free_diag(m); }

ld vector_get_elemB( vector *v, ui i ) { return elemB( v->size, i ); }
void vector_print_col( vector *v );

vector *vector_subtract( vector *a, vector *b );
ld vector_norm_inf( vector *v );

vector *vector_new_full( ui n );
vector *vector_new_full_copy( vector *v );
ld vector_get_full( vector *v, ui i );
bool vector_put_full( vector *v, ui i, ld a );
void vector_free_full( vector *v );

int main( int argc, const char * argv[] )
{
/*	// example from http://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method#An_example_for_the_matrix_version
	matrix *A = matrix_new_full( 2 );
	A->put( A,1,1, 16 );
	A->put( A,1,2, 3 );
	A->put( A,2,1, 7 );
	A->put( A,2,2, -11 );
	vector *b = vector_new_full( 2 );
	b->put( b,1, 11 );
	b->put( b,2, 13 );
	vector *z = vector_new_full( 2 );
	z->put( z,1, 1 );
	z->put( z,2, 1 );
	printf( "A=\n" );
	matrix_print(A);
	printf( "b=\n" );
	vector_print_col( b );
	struct gauss_info g = solve_gauss( A, b, z, 1e-8L, 20 );
	vector *x_gauss = g.x;
	vector *b_gauss = g.b;
	vector *r_gauss = g.r;
	ld n_gauss = g.n;
	printf( "gauss: solution x=\n" );
	vector_print_col( x_gauss );
	printf( "gauss: recalculated b=\n" );
	vector_print_col( b_gauss );
	printf( "gauss: residual r=\n" );
	vector_print_col( r_gauss );
	printf( "infinity norm of gauss residual: n = % 8.6Le  = %43.40Lf\n", n_gauss, n_gauss );
	printf( "gauss iterations: %u\n", g.i );
	return 0;
*/	//ui n, d=1;
	ui n;
	matrix A={ &matrix_get_elemA, &matrix_put_fail, n, NULL };
	vector b={ &vector_get_elemB, &vector_put_fail, n, NULL };
	vector z={ &vector_get_zero, &vector_put_fail, n, NULL };
	//for( n = 1; n <= 1000000; n += d ) {
	for( n = 10; n <= 10000; n *= 10 ) {
		A.size = n;
		b.size = n;
		z.size = n;
		printf( "\nn=%u\n", A.size );
		//printf( "\u2016A\u2016_\u221E = %43.40Lf\n", norm_inf(&A) );
		//for( i = 1; i <= n; i++ ) { printf( "n=%u, i=%u, b_i = %43.40Lf\n", n, i, elemB(n,i) ); }
		if ( n <= 0 ) {
			printf( "A=\n" );
			matrix_print(&A);
			printf( "b=\n" );
			vector_print_col( &b );
		}
		vector *x_crout = solve_crout( &A, &b );
		vector *b_crout = matrix_apply_vector_col( &A, x_crout );
		vector *r_crout = vector_subtract( &b, b_crout );
		ld n_crout = vector_norm_inf( r_crout );
		if ( n <= 0 ) {
			//printf( "L=\n" );
			//matrix_print( LU[0] );
			//printf( "U=\n" );
			//matrix_print( LU[1] );
			//printf( "z=\n" );
			//vector_print_col( z );
			printf( "crout: solution x=\n" );
			vector_print_col( x_crout );
			printf( "crout: recalculated b=\n" );
			vector_print_col( b_crout );
			printf( "crout: residual r=\n" );
			vector_print_col( r_crout );
		}
		printf( "infinity norm of crout residual: n = % 8.6Le  = %43.40Lf  (not iterative)\n", n_crout, n_crout );
		vector_free_full( x_crout );
		vector_free_full( b_crout );
		vector_free_full( r_crout );
		//struct gauss_info g = solve_gauss( &A, &b, &z, r_crout*2, 1000000 );
		struct gauss_info_fast f = solve_gauss_fast( &A, &b, &z, 1e-8L, 10000000 );
		printf( "infinity norm of gauss residual: n = % 8.6Le  = %43.40Lf  (%u iterations)\n", f.n, f.n, f.i );
/*		struct gauss_info g = solve_gauss( &A, &b, &z, 1e-8L, 1000000 );
		vector *x_gauss = g.x;
		vector *b_gauss = g.b;
		vector *r_gauss = g.r;
		ld n_gauss = g.n;
		if ( n <= 0 ) {
			printf( "gauss: solution x=\n" );
			vector_print_col( x_gauss );
			printf( "gauss: recalculated b=\n" );
			vector_print_col( b_gauss );
			printf( "gauss: residual r=\n" );
			vector_print_col( r_gauss );
		}
		//		printf( "infinity norm of crout residual: n = % 8.6Le  = %43.40Lf  (not iterative)\n", n_crout, n_crout );
		printf( "infinity norm of gauss residual: n = % 8.6Le  = %43.40Lf  (%u iterations)\n", n_gauss, n_gauss, g.i );
		vector_free_full( x_crout );
		vector_free_full( b_crout );
		vector_free_full( r_crout );
		vector_free_full( x_gauss );
		vector_free_full( b_gauss );
		vector_free_full( r_gauss );
*/		//if( d <= n/10 ) { d *= 10; }
	}
	printf( "\n" );
	printf( "1. Crout factorization has much smaller error for each n.\n" );
	printf( "2. As n increases, Crout factorization error increases very slowly,\n" );
	printf( "               but Gauss-Seidel iteration error gets to the same target\n" );
	printf( "3. As n increases, Gauss-Seidel iterations required increases very fast\n" );
	printf( "       - nearly 2 orders of magnitude for every order of magnitude of n.\n" );
	printf( "\n" );
	printf( "This is using the \"fast\" implementation of Gauss-Seidel, which is specific to this matrix A.\n" );
	printf( "On my computer, at n=1000, the slow (and more general) implementation takes over 5 hours,\n" );
	printf( "while the fast implementation takes about 37 seconds.\n" );
	printf( "(The slow implementation is still in the code.)\n" );
	printf( "\n" );
    return 0;
}

ld elemA( ui n, ui i, ui j ) {
	ld r = 0;
	if( i == j ) { r = 2; }
	else if( udiff(i,j) == 1 ) { r = -1; }
	return r;
}

ld elemB( ui n, ui i ) {
	ld r = ((ld)i)/powl(((ld)n)+1,(ld)4);
	if( 1 == i ) { r += 1; }
	else if( n == i ) { r += 6; }
	return r;
}

vector *solve_crout( matrix *A, vector *b ) {
	if( A->size != b->size ) { crash(); } // fail early, fail loud
	matrix **LU = factor_crout( A, 1 );
	vector *z = fwdsub( LU[0], b );
	vector *x = baksub( LU[1], z );
	matrix_free_lower(LU[0]);
	matrix_free_upper(LU[1]);
	free(LU);
	vector_free_full( z );
	return x;
}

matrix **factor_crout( matrix *A, ui k ) {
	matrix **r = malloc( 3 * sizeof(matrix *) );
	matrix *L = r[0] = matrix_new_lower( A->size, k );
	matrix *U = r[1] = matrix_new_upper( A->size, k );
	r[2] = NULL;
	ui n=A->size, s, t, v;
	ld c, d;
	// modeled after http://mathfaculty.fullerton.edu/mathews/n2003/CholeskyMod.html
	for( s = 1; s <= n; s++ ) { // step along diagonal
		for( t = s; t <= n; t++ ) { // step away from diagonal*/
			c = d = 0;
			for( v = umax(1,usub(s,k)); v < s ; v++ ) {
				c += L->get(L,t,v)*U->get(U,v,s);
				d += L->get(L,s,v)*U->get(U,v,t);
			}
			//printf( "crout s=%u, t=%u; c=% 8.6Le, d=% 8.6Le\n", s, t, c, d );
			L->put( L,t,s, (A->get(A,t,s) - c) );
			U->put( U,s,t, (A->get(A,s,t) - d) / L->get(L,s,s) );
		}
	}
	return r;
}

vector *fwdsub( matrix *L, vector *b ) {
	if( L->size != b->size ) { crash(); } // fail early, fail loud
	ui n = L->size, i, j;
	ld c;
	vector *r = vector_new_full( n );
	for( i = 1; i <= n; i++ ) {
		c = 0;
		for( j = 1; j < i; j++ ) {
			c += L->get(L,i,j) * r->get(r,j);
		}
		r->put( r,i, (b->get(b,i) - c) / L->get(L,i,i) );
	}
	return r;
}

vector *baksub( matrix *U, vector *z ) {
	if( U->size != z->size ) { crash(); } // fail early, fail loud
	ui n = U->size, i, j;
	ld c;
	vector *r = vector_new_full( n );
	for( i = n; i >= 1; i-- ) {
		c = 0;
		for( j = n; j > i; j-- ) {
			c += U->get(U,i,j) * r->get(r,j);
		}
		r->put( r,i, (z->get(z,i) - c) / U->get(U,i,i) );
	}
	return r;
}

struct gauss_info_fast solve_gauss_fast( matrix *A, vector *b, vector *x0, ld tolerance, ui ceiling ) {
	if( A->size != b->size || A->size != x0->size ) { crash(); } // fail early, fail loud
	if( A->get != matrix_get_elemA ) { crash(); } // fast is specialized
	if( A->size < 3 ) { crash(); } // fast is specialized
	ui n=A->size, i, j;
	ld nn, np, s;
	ld *xn, *xp, *xt, *bc;
	xn = malloc( n * sizeof(ld) ); //next x (calculated in main loop)
	xp = malloc( n * sizeof(ld) ); //previous x (copied from x0)
	bc = malloc( n * sizeof(ld) ); //copied b
	for( i = 0; i < n; i++ ) {
		xp[i] = x0->get(x0,i+1);
		bc[i] = b->get(b,i+1);
	}
	nn = np = 10e10L*tolerance;
	i = 0;
	while( i < ceiling && nn > tolerance && nn <= np ) {
		xn[0] = (bc[0] + xp[1]) / 2;                    // new x element 0
		for( j = 1; j < n-1; j++ ) {
			xn[j] = (bc[j] + (xn[j-1] + xp[j+1]) ) / 2; // new x element j
		}
		xn[j] = (bc[n-1] + xn[n-2]) / 2;                // new x element n
		np = nn; nn = 0; // save previous norm, reset new norm
		s = 2 * xn[0] - xn[1];               // estimated b element 0
		s = bc[0]>s ? bc[0]-s : s-bc[0];     // absolute residual element 0
		if( s > nn ) { nn = s; }             // norm of elements 0 through 0
		for( j = 1; j < n-1; j++ ) {
			s = 2*xn[j] - xn[j-1] - xn[j+1]; // estimated b element j
			s = bc[j]>s ? bc[j]-s : s-bc[j]; // absolute residual element j
			if( s > nn ) { nn = s; }         // norm of elements 0 through j
		}
		j = n-1;
		s = -xn[j-1] + 2 * xn[j];            // estimated b element n
		s = bc[j]>s ? bc[j]-s : s-bc[j];     // absolute residual element n
		if( s > nn ) { nn = s; }             // norm of elements 0 through n
		++i; xt = xp; xp = xn; xn = xt; //cycle
	}
	struct gauss_info_fast g = { nn, i };
	return g;
}

struct gauss_info solve_gauss( matrix *A, vector *b, vector *x0, ld tolerance, ui ceiling ) {
	if( A->size != b->size || A->size != x0->size ) { crash(); } // fail early, fail loud
	ui n=A->size, i=0;
	ld m = 10e10L*tolerance;
	vector *t = vector_new_full( n ); // because we need two writable vectors
	struct gauss_info g;
	g.x = vector_new_full_copy( x0 );
	g.n = m;
//	printf( "gauss iteration %u: %p (%p)\n", i, g.x, t ); vector_print_col( g.x );
	while( i < ceiling && g.n > tolerance && g.n <= m ) {
		x0 = g.x; // keep a pointer to the previous iteration
		m = g.n; // previous norm should be bigger then new norm
		g.x = gauss_iterate( t, A, b, g.x ); // updates t using g.x as previous, then reassigns g.x
		g.b = matrix_apply_vector_col( A, g.x );
		g.r = vector_subtract( b, g.b );
		g.n = vector_norm_inf( g.r );
		g.i = ++i;
		t = x0; // previous becomes next
//		printf( "gauss iteration %u: %p (%p)\n", g.i, g.x, t ); vector_print_col( g.x );
	}
	return g;
}

// updates xn, uses xp for previous
vector *gauss_iterate( vector *xn, matrix *A, vector *b, vector *xp ) {
	ui n = A->size, i, j;
	ld d;
	for( i = 1; i <= n; i++ ) {
		d = 0;
		for( j = 1; j < i; j++ ) {
			d += A->get(A,i,j) * xn->get(xn,j);
		}
		// j=i is skipped
		for( j = i+1; j <= n; j++ ) {
			d += A->get(A,i,j) * xp->get(xp,j);
		}
		xn->put( xn,i, (b->get(b,i) - d) / A->get(A,i,i) );
	}
	return xn;
}


void matrix_print( matrix *m ) {
	ui i=1, j=1;
	switch( m->size ) {
		case 0:
			break;
		case 1:
			printf( "\uFF3B " );
			for(; j <= m->size; j++) { printf( " % 8.6Le ", m->get(m,i,j) ); }
			printf( " \uFF3D\n" );
			break;
		default:
			printf( "\u23A1 " );
			for(; j <= m->size; j++) { printf( " % 8.6Le ", m->get(m,i,j) ); }
			printf( " \u23A4\n" );
			for( i = 2; i < m->size; i++ ) {
				printf( "\u23A2 " );
				for( j = 1; j <= m->size; j++) { printf( " % 8.6Le ", m->get(m,i,j) ); }
				printf( " \u23A5\n" );
			}
			printf( "\u23A3 " );
			for( j = 1; j <= m->size; j++) { printf( " % 8.6Le ", m->get(m,i,j) ); }
			printf( " \u23A6\n" );
	}
}

vector *matrix_apply_vector_col( matrix *A, vector *x ) {
	if( A->size != x->size ) { crash(); } // fail early, fail loud
	ui n = A->size, i, j;
	ld b;
	vector *r = vector_new_full( n );
	for( i = 1; i <= n; i++ ) {
		b = 0;
		for( j = 1; j <= n; j++ ) {
			b += A->get(A,i,j) * x->get(x,j);
		}
		r->put( r,i, b );
	}
	return r;
}

ld matrix_norm_inf( matrix *A ) {
	ld r = 0, s = 0;
	ui n = A->size, i = 1, j = 1;
	for(; i <= n; i++ ) {
		for(; j <= n; j++ ) {
			s += A->get( A, i, j );
			//printf( "n=%u, i=%u, j=%u, s_(i,j) = %43.40Lf\n", n, i, j, s );
		}
		s = fabsl(s);
		//printf( "n=%u, i=%u, \u2211_i = %43.40Lf\n", n, i, s );
		if( s > r ) { r = s; }
		s = 0;
	}
	return r;
}

matrix *matrix_new_full( ui n ) {
	matrix *r = malloc( sizeof(matrix) );
	r->get = matrix_get_full;
	r->put = matrix_put_full;
	r->size = n;
	r->data = calloc( n*n, sizeof(ld) );
	return r;
}

ld matrix_get_full( matrix *m, ui i, ui j ) {
	ui n = m->size, f = n*(i-1) + j;
	if( f > n*n ) { crash(); }
	return ((ld*)m->data)[f-1];
}

bool matrix_put_full( matrix *m, ui i, ui j, ld a ) {
	ui n = m->size, f = n*(i-1) + j;
	if( f > n*n ) { crash(); }
	((ld*)m->data)[f-1] = a;
	return true;
}

void matrix_free_full( matrix *m ) {
	free( m->data );
	free( m );
}

matrix *matrix_new_diag( ui n, ui k ) { // n-by-n with diagional and k sub/superdiagonals
	ui s = (k+1)*(2*n-k)/2; // number of elements to store; n + (n-1) + ... + (n-k)
	ld *e = calloc( s, sizeof(ld) );
	struct matrix_data_diag *d = malloc( sizeof(struct matrix_data_diag) );
	matrix *m = malloc( sizeof(matrix) );
	d->k = umin( k, n-1 );
	d->s = s; // just in case
	d->e = e;
	m->get = &matrix_get_zero;
	m->put = &matrix_put_fail;
	m->size = n;
	m->data = d;
	return m;
}

void matrix_free_diag( matrix *m ) {
	free( ((struct matrix_data_diag *)m->data)->e );
	free( m->data );
	free( m );
}

matrix *matrix_new_lower( ui n, ui k ) {
	matrix *r = matrix_new_diag( n, k );
	r->get = &matrix_get_lower;
	r->put = &matrix_put_lower;
	return r;
}

ld matrix_get_lower( matrix *m, ui i, ui j ) {
	if(0==i){crash();}
	ld r=0;
	if( j <= i ) { // j > i is upper, keep r=0
		struct matrix_data_diag *d = m->data;
		ui k = i-j, l=j; // lth element in kth subdiagonal
		if( k <= d->k ) { // if subdiagonal is too low, keep r=0
			// if the element array were square, each k would start at n*k
			// to avoid unused elements, it's n*k - 0 - 1 - 2 - ... - (k-1)
			ui f = m->size*k - k*(k-1)/2 + l;
//			printf( "get_lower( n=%u, i=%u, j=%u ): element %u\n", m->size, i, j, f );
			if( f > d->s ) { crash(); } // fail early, fail loud
			r = d->e[usub(f,1)]; // -1 because e starts at 0
		}
//		else{ printf( "get_lower( n=%u, i=%u, j=%u ): k=%u is too low\n", m->size, i, j, k ); }
	}
//	else { printf( "get_lower( n=%u, i=%u, j=%u ): upper\n", m->size, i, j ); }
	return r;
}

bool matrix_put_lower( matrix *m, ui i, ui j, ld a ) {
	bool r=false;
	if( j <= i ) { // j > i is upper, keep r=false
		struct matrix_data_diag *d = m->data;
		ui k = i-j, l=j; // lth element in kth subdiagonal
		if( k <= d->k ) { // if subdiagonal is too low, keep r=false
			// if the element array were square, each k would start at n*k
			// to avoid unused elements, it's n*k - 0 - 1 - 2 - ... - (k-1)
			ui f = m->size*k - k*(k-1)/2 + l;
//			printf( "put_lower( n=%u, i=%u, j=%u, % 8.6Le ): element %u\n", m->size, i, j, a, f );
			if( f > d->s ) { crash(); } // fail early, fail loud
			if( 0 != d->e[usub(f,1)] ) { crash(); } // usually wrong
			d->e[usub(f,1)] = a; // -1 because e starts at 0
			r = true;
		}
//		else{ printf( "put_lower( n=%u, i=%u, j=%u, % 8.6Le ): k=%u is too low\n", m->size, i, j, a, k ); }
	}
//	else { printf( "put_lower( n=%u, i=%u, j=%u, % 8.6Le ): upper\n", m->size, i, j, a ); }
	return r;
}

matrix *matrix_new_upper( ui n, ui k ) {
	matrix *r = matrix_new_diag( n, k );
	r->get = &matrix_get_upper;
	r->put = &matrix_put_upper;
	return r;
}

ld matrix_get_upper( matrix *m, ui i, ui j ) {
	if(0==i){crash();}
	ld r=0;
	if( j >= i ) { // j < i is lower, keep r=0
		struct matrix_data_diag *d = m->data;
		ui k = j-i, l=i; // lth element in kth superdiagonal
		if( k <= d->k ) { // if superdiagonal is too high, keep r=0
			// if the element array were square, each k would start at n*k
			// to avoid unused elements, it's n*k - 0 - 1 - 2 - ... - (k-1)
			ui f = m->size*k - k*(k-1)/2 + l;
//			printf( "get_upper( n=%u, i=%u, j=%u ): element %u\n", m->size, i, j, f );
			if( f > d->s ) { crash(); } // fail early, fail loud
			r = d->e[usub(f,1)]; // -1 because e starts at 0
		}
//		else{ printf( "get_upper( n=%u, i=%u, j=%u ): k=%u is too high\n", m->size, i, j, k ); }
	}
//	else { printf( "get_upper( n=%u, i=%u, j=%u ): lower\n", m->size, i, j ); }
	return r;
}

bool matrix_put_upper( matrix *m, ui i, ui j, ld a ) {
	bool r=false;
	if( j >= i ) { // j < i is lower, keep r=false
		struct matrix_data_diag *d = m->data;
		ui k = j-i, l=i; // lth element in kth superdiagonal
		if( k <= d->k ) { // if superdiagonal is too high, keep r=false
			// if the element array were square, each k would start at n*k
			// to avoid unused elements, it's n*k - 0 - 1 - 2 - ... - (k-1)
			ui f = m->size*k - k*(k-1)/2 + l;
//			printf( "put_upper( n=%u, i=%u, j=%u, % 8.6Le ): element %u\n", m->size, i, j, a, f );
			if( f > d->s ) { crash(); } // fail early, fail loud
			if( 0 != d->e[usub(f,1)] ) { crash(); } // usually wrong
			d->e[usub(f,1)] = a; // -1 because e starts at 0
			r = true;
		}
//		else{ printf( "put_upper( n=%u, i=%u, j=%u, % 8.6Le ): k=%u is too high\n", m->size, i, j, a, k ); }
	}
//	else { printf( "put_upper( n=%u, i=%u, j=%u, % 8.6Le ): lower\n", m->size, i, j, a ); }

	return r;
}

void vector_print_col( vector *v ) {
	ui i=1;
	switch( v->size ) {
		case 0:
			break;
		case 1:
			printf( "\uFF3B  % 8.6Le  \uFF3D\n", v->get(v,i) );
			break;
		default:
			printf( "\u23A1  % 8.6Le  \u23A4\n", v->get(v,i) );
			for( i = 2; i < v->size; i++ ) {
				printf( "\u23A2  % 8.6Le  \u23A5\n", v->get(v,i) );
			}
			printf( "\u23A3  % 8.6Le  \u23A6\n", v->get(v,i) );
	}
}

vector *vector_subtract( vector *a, vector *b ) {
	if( a->size != b->size ) { crash(); } // fail early, fail loud
	ui n = a->size, i;
	vector *r = vector_new_full( n );
	for( i = 1; i <= n; i++ ) {
		r->put( r,i, a->get(a,i) - b->get(b,i) );
	}
	return r;
}

ld vector_norm_inf( vector *v ) {
	ld r = 0, s;
	ui n = v->size, i;
	for( i = 1; i <= n; i++ ) {
		s = fabsl( v->get( v, i ) );
		if( s > r ) { r = s; }
	}
	return r;
}

vector *vector_new_full( ui n ) {
	vector *r = malloc( sizeof(vector) );
	r->get = &vector_get_full;
	r->put = &vector_put_full;
	r->size = n;
	r->data = calloc( n, sizeof(ld) );
	return r;
}

vector *vector_new_full_copy( vector *v ) {
	ui n = v->size, i;
	vector *r = malloc( sizeof(vector) );
	r->get = &vector_get_full;
	r->put = &vector_put_full;
	r->size = n;
	r->data = malloc( n * sizeof(ld) );
	for( i = 1; i <= n; i++ ) {
		r->put( r,i, v->get(v,i) );
	}
	return r;
}

ld vector_get_full( vector *v, ui i ) {
	ld r = 0;
	if( i != 0 && i <= v->size ) { r = ((ld*)v->data)[i-1]; }
	return r;
}

bool vector_put_full( vector *v, ui i, ld a ) {
	bool r = false;
	if( i != 0 && i <= v->size ) {
		((ld*)v->data)[i-1] = a;
		r = true;
	}
	return r;
}

void vector_free_full( vector *v ) {
	free( v->data);
	free( v );
}





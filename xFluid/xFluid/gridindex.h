#pragma once


struct GridIndex3d {
public:

	GridIndex3d() : i(0), j(0), k(0) {}
	GridIndex3d(int ii, int jj, int kk) : i(ii), j(jj), k(kk) {}

	bool operator==(const GridIndex3d &other) const {
		return i == other.i &&
			j == other.j &&
			k == other.k;
	}

	bool operator!=(const GridIndex3d &other) const {
		return i != other.i ||
			j != other.j ||
			k != other.k;
	}

public:
	int i, j, k;
};


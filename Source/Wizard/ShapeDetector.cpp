// Fill out your copyright notice in the Description page of Project Settings.

#include "ShapeDetector.h"
#include <Eigen/Dense>



FString UShapeDetector::DetectShape()
{
	return "Hey it's triangle";
}

float UShapeDetector::getCircle2fResemblance(const TArray<FVector2f>& points) {
	FVector2f center = FVector2f::ZeroVector;

	for (const FVector2f& p : points) {
		center += p;
	}
	center /= points.Num();

	//UE_LOG(LogTemp, Display, TEXT("Calculated center: %f %f"), center.X, center.Y);

	float avgDist = 0;
	TArray<float> distances;
	for (const FVector2f& p : points) {
		float dist = FVector2f::Distance(p, center);
		distances.Add(dist);
		avgDist += dist;
	}
	avgDist /= points.Num();

	//UE_LOG(LogTemp, Display, TEXT("AVG dist is: %f"), avgDist	);


	float avgDeviation = 0;
	for (int i=0; i<points.Num(); i++) {
		avgDeviation += abs(distances[i] - avgDist);
	}
	avgDeviation /= points.Num();

	//UE_LOG(LogTemp, Display, TEXT("AVG Deviation is: %f"), avgDeviation);

	
	return 1 - avgDeviation / avgDist;
}

bool UShapeDetector::IsCircle2f(const TArray<FVector2f>& points, float precision) {
	FVector2f center = FVector2f::ZeroVector;

	for (const FVector2f& p : points) {
		center += p;
	}
	center /= points.Num();

	//UE_LOG(LogTemp, Display, TEXT("Calculated center: %f %f"), center.X, center.Y);

	float avgDist = 0;
	TArray<float> distances;
	for (const FVector2f& p : points) {
		float dist = FVector2f::Distance(p, center);
		distances.Add(dist);
		avgDist += dist;
	}
	avgDist /= points.Num();

	//UE_LOG(LogTemp, Display, TEXT("AVG dist is: %f"), avgDist	);


	float avgDeviation = 0;
	for (int i=0; i<points.Num(); i++) {
		avgDeviation += abs(distances[i] - avgDist);
	}
	avgDeviation /= points.Num();

	//UE_LOG(LogTemp, Display, TEXT("AVG Deviation is: %f"), avgDeviation);

	
	return getCircle2fResemblance(points) + precision > 1;
}


FPlane UShapeDetector::ComputeBestFitPlane(const TArray<FVector3f>& points) {
	if (points.Num() < 3) {
		return FPlane(); 
	}

	FVector3f centroid = FVector3f::ZeroVector;
	for (const FVector3f& p : points) {
		centroid += p;
	}
	centroid /= points.Num();

	Eigen::Matrix3d covarianceMatrix;
	covarianceMatrix.setZero();

	for (const FVector3f& p : points) {
		FVector3f d = p - centroid;
		covarianceMatrix(0, 0) += d.X * d.X;
		covarianceMatrix(0, 1) += d.X * d.Y;
		covarianceMatrix(0, 2) += d.X * d.Z;
		covarianceMatrix(1, 1) += d.Y * d.Y;
		covarianceMatrix(1, 2) += d.Y * d.Z;
		covarianceMatrix(2, 2) += d.Z * d.Z;
	}
	covarianceMatrix(1, 0) = covarianceMatrix(0, 1);
	covarianceMatrix(2, 0) = covarianceMatrix(0, 2);
	covarianceMatrix(2, 1) = covarianceMatrix(1, 2);

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covarianceMatrix);
	Eigen::Vector3d normalEigenVector = solver.eigenvectors().col(0); // Собственный вектор с наименьшим значением

	FVector3f normal(normalEigenVector(0), normalEigenVector(1), normalEigenVector(2));
	normal.Normalize();

	float w = -(normal.X * centroid.X + normal.Y * centroid.Y + normal.Z * centroid.Z);

		
	return FPlane(normal.X, normal.Y, normal.Z, w); // w - do i need it?
}
								  
TArray<FVector2f> UShapeDetector::ProjectPointsToPlane(const TArray<FVector3f>& points, const FPlane& plane) {
	TArray<FVector2f> projectedPoints;

	// Проверяем нормаль плоскости
	FVector3f normal(plane.X, plane.Y, plane.Z);
	normal.Normalize(); 

	// Выбираем оси для 2D координат (основание плоскости)
	FVector3f basisX, basisY;
	FVector3f upVector(0, 0, 1); // Стандартный up-вектор
	
	if (FMath::Abs(normal.Z) > 0.9f) {
		upVector = FVector3f(1, 0, 0); // Если нормаль почти вертикальна, меняем up-вектор
	}

	basisX = FVector3f::CrossProduct(upVector, normal).GetSafeNormal();
	basisY = FVector3f::CrossProduct(normal, basisX).GetSafeNormal();

	// Проецируем каждую точку
	for (const FVector3f& point : points) {
		// Находим проекцию точки на плоскость
		float distance = FVector3f::DotProduct(point, normal) - plane.W;
		FVector3f projectedPoint = point - normal * distance;

		// Преобразуем в 2D координаты
		float u = FVector3f::DotProduct(projectedPoint, basisX);
		float v = FVector3f::DotProduct(projectedPoint, basisY);

		projectedPoints.Add(FVector2f(u, v));
	}

	return projectedPoints;
}

bool UShapeDetector::IsPointNearLine(const FVector2f& p, const FVector2f& p1, const FVector2f& p2, float epsilon) {
	if (p1 == p2) return false; // Точки совпадают — некорректная линия

	// Вектор направления линии
	FVector2f dir = p2 - p1;
	dir.Normalize();

	// Вектор из p1 в точку p
	FVector2f v = p - p1;

	// Проекция v на направление линии
	float projection = FVector2f::DotProduct(v, dir);

	// Ближайшая точка на линии
	FVector2f closestPoint = p1 + dir * projection;

	// Проверяем расстояние до линии
	return FVector2f::DistSquared(p, closestPoint) < epsilon * epsilon;
}


// Main API
bool UShapeDetector::IsCircle3f(const TArray<FVector3f>& points, float precision) {
	FPlane plane = ComputeBestFitPlane(points);
	TArray<FVector2f> points2f= ProjectPointsToPlane(points, plane);

	return IsCircle2f(points2f, precision);;
}

float UShapeDetector::getCircle3fResemblance(const TArray<FVector3f>& points) {
	FPlane plane = ComputeBestFitPlane(points);
	TArray<FVector2f> points2f= ProjectPointsToPlane(points, plane);

	return getCircle2fResemblance(points2f);
}

bool UShapeDetector::IsTriangle(const TArray<FVector2f>& points, float epsilon, int maxIterations) {
	if (points.Num() < 3) return false;

	TArray<FVector2f> remainingPoints = points;
	TArray<TPair<FVector2f, FVector2f>> lines;

	for (int i = 0; i < 3*maxIterations; i++) { // Ищем 3 линии
		if (remainingPoints.Num() < 2) return false;

		// RANSAC: выбираем 2 случайные точки
		int index1 = FMath::RandRange(0, remainingPoints.Num() - 1);
		int index2;
		do {
			index2 = FMath::RandRange(0, remainingPoints.Num() - 1);
		} while (index1 == index2);

		FVector2f p1 = remainingPoints[index1];
		FVector2f p2 = remainingPoints[index2];

		// Проверяем, сколько точек попадает на линию
		TArray<FVector2f> inliers;
		for (const FVector2f& p : remainingPoints) {
			if (IsPointNearLine(p, p1, p2, epsilon)) {
				inliers.Add(p);
			}
		}

		// Если нашли достаточно (хотя бы 20%) точек, фиксируем линию
		if (inliers.Num() >= points.Num()/5) {
			lines.Add(TPair<FVector2f, FVector2f>(p1, p2));

			// Удаляем использованные точки
			for (const FVector2f& p : inliers) {
				remainingPoints.Remove(p);
			}
		}
	}

	// Если нашли 3 линии, то это треугольник
	return lines.Num() == 3;
}

/*TArray<FVector2f> UShapeDetector::GetPointsCommonLines(const TArray<FVector2f>& points, float epsilon, int maxIterations) {
	//if (points.Num() < 3) return false;

	TArray<FVector2f> remainingPoints = points;
	TArray<FVector2f> lines;

	for (int i = 0; i < maxIterations; i++) { // Ищем 3 линии
		//if (remainingPoints.Num() < 2) return false;

		// RANSAC: выбираем 2 случайные точки
		int index1 = FMath::RandRange(0, remainingPoints.Num() - 1);
		int index2;
		do {
			index2 = FMath::RandRange(0, remainingPoints.Num() - 1);
		} while (index1 == index2);

		FVector2f p1 = remainingPoints[index1];
		FVector2f p2 = remainingPoints[index2];

		// Проверяем, сколько точек попадает на линию
		TArray<FVector2f> inliers;
		for (const FVector2f& p : remainingPoints) {
			if (IsPointNearLine(p, p1, p2, epsilon)) {
				inliers.Add(p);
			}
		}

		// Если нашли достаточно (хотя бы 20%) точек, фиксируем линию
		if (inliers.Num() >= points.Num()/5) {
			lines.Add(FVector2f(p1 - p2));

			// Удаляем использованные точки
			for (const FVector2f& p : inliers) {
				remainingPoints.Remove(p);
			}
		}
	}

	return lines;
}*/

TArray<FVector2f> UShapeDetector::GetPointsCommonLines(const TArray<FVector2f>& points, float epsilon, int maxIterations) {
	// Копируем все точки
	TArray<FVector2f> remainingPoints = points;
	TArray<FVector2f> lines;

	// Пытаемся найти 3 линии
	int pointsProcessed = 0;
	for (int i = 0; i < maxIterations; i++) {
		// Если слишком мало точек, выходим
		if (remainingPoints.Num() < 2) break;

		// RANSAC: выбираем 2 случайные точки
		int index1 = FMath::RandRange(0, remainingPoints.Num() - 1);
		int index2;
		do {
			index2 = FMath::RandRange(0, remainingPoints.Num() - 1);
		} while (index1 == index2);

		FVector2f p1 = remainingPoints[index1];
		FVector2f p2 = remainingPoints[index2];

		// Проверяем, сколько точек попадает на линию
		TArray<FVector2f> inliers;
		for (const FVector2f& p : remainingPoints) {
			if (IsPointNearLine(p, p1, p2, epsilon)) {
				inliers.Add(p);
			}
		}

		// Если нашли достаточно (хотя бы 16.6%) точек, фиксируем линию
		if (inliers.Num() >= points.Num() / 6) {
			// Добавляем линию (вектор)
			lines.Add(p2 - p1);
			pointsProcessed += inliers.Num(); 

			// Удаляем использованные точки
			for (const FVector2f& p : inliers) {
				remainingPoints.Remove(p);
			}
		}

		if (pointsProcessed * 0.8 > points.Num()) {
			//UE_LOG(LogTemp, Display, TEXT("Points processed break at: %i"), i);

			break;
		}
	}

	return lines;
}

TArray<FVector3f> UShapeDetector::FilterPointsFromSameOrigin(const TArray<FVector3f>& points, float radius, int count) {
	TArray<FVector3f> filteredPoints;
	float radiusSquared = radius * radius; 

	for (const FVector3f& point : points) {
		int neighbors = 0;

		for (const FVector3f& otherPoint : filteredPoints) {
			if (&point != &otherPoint && FVector3f::DistSquared(point, otherPoint) <= radiusSquared) {
				neighbors++;
				if (neighbors >= count) break; 
			}
		}

		if (neighbors < count) {
			filteredPoints.Add(point);
		}
	}

	return filteredPoints;
}

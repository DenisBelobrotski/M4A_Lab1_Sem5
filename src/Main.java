public class Main {
    public static void main(String[] args) {
        double intervalBottom = 0;
        double intervalUpper = 1;
        double accuracy = Math.pow(10, -5);
        int nodeNumGauss = 4;
        System.out.println("Формула левых прямоугольников по правилу Рунге:");
        System.out.println("Значение интеграла: " + calcIntegralRungeLeftRectangles(intervalBottom, intervalUpper, accuracy));
        System.out.println();
        System.out.println("Формула средних прямоугольников:");
        System.out.println("Значение интеграла: " + calcIntegralMiddleRectangles(intervalBottom, intervalUpper, accuracy));
        System.out.println();
        System.out.println("Формула Симпсона:");
        System.out.println("Значение интеграла: " + calcIntegralSimpson(intervalBottom, intervalUpper, accuracy));
        System.out.println();
        System.out.println("Квадратурная формула типа Гаусса:");
        System.out.println("Значение интеграла: " + calcIntegralGauss(intervalBottom, intervalUpper, nodeNumGauss));
    }

    private static double calcIntegralRungeLeftRectangles(double intervalBottom, double intervalUpper, double accuracy) {
        double prevStep;
        double curStep = intervalUpper - intervalBottom;
        double prevIntegral;
        double curIntegral = calcIntegralLeftRectangles(intervalBottom, intervalUpper, curStep);
        double rungeCoef;
        do {
            prevStep = curStep;
            curStep /= 2;
            prevIntegral = curIntegral;
            curIntegral = calcIntegralLeftRectangles(intervalBottom, intervalUpper, curStep);
            rungeCoef = calcRungeResidual(prevStep, curStep, prevIntegral, curIntegral, 1);
        } while (Math.abs(rungeCoef) > accuracy);
        System.out.println("Шаг: " + curStep);
        System.out.println("Точность: " + accuracy);
        System.out.println("Количество узлов: " + (new Double(Math.ceil((intervalUpper - intervalBottom) / curStep))).intValue());
        return curIntegral + (curIntegral - prevIntegral) / (Math.pow(prevStep / curStep, 1) - 1);
    }

    private static double calcRungeResidual(double prevStep, double curStep, double prevIntegral, double curIntegral, double power) {
        return (curIntegral - prevIntegral) / (1 - Math.pow(curStep / prevStep, power));
    }

    private static double calcIntegralLeftRectangles(double intervalBottom, double intervalUpper, double step) {
        Double doubleNodeNum = (intervalUpper - intervalBottom) / step;
        int nodeNum = doubleNodeNum.intValue();
        double sum = 0;
        for (int i = 0; i < nodeNum - 1; i++) {
            sum += calcFunction(intervalBottom + i * step);
        }
        return step * sum;
    }

    private static double calcFunction(double x) {
        return 1 / (1 + Math.pow(x, 3));
    }

    private static double calcIntegralMiddleRectangles(double intervalBottom, double intervalUpper, double accuracy) {
        int nodeNum = calcNodeNumMiddleRectangle(intervalBottom, intervalUpper, accuracy);
        double step = (intervalUpper - intervalBottom) / nodeNum;
        double sum = 0;
        System.out.println("Количество узлов: " + nodeNum);
        System.out.println("Шаг: " + step);
        for (int i = 0; i <= nodeNum - 1; i++) {
            sum += calcFunction(intervalBottom + (i + 0.5) * step);
        }
        return step * sum;
    }

    private static int calcNodeNumMiddleRectangle(double intervalBottom, double intervalUpper, double accuracy) {
        Double result = Math.ceil(Math.sqrt(Math.pow(intervalUpper - intervalBottom, 3) * calcSecondDerivative(intervalUpper) / (24 * accuracy)));
        return result.intValue();
    }

    private static double calcSecondDerivative(double x) {
        return 6 * x * (2 * Math.pow(x, 3) - 1) / Math.pow(Math.pow(x, 3) + 1, 3);
    }

    private static double calcIntegralSimpson(double intervalBottom, double intervalUpper, double accuracy) {
        int nodeNum = calcNodeNumSimpson(intervalBottom, intervalUpper, accuracy);
        double step = (intervalUpper - intervalBottom) / nodeNum;
        double sum1 = 0;
        double sum2 = 0;
        System.out.println("Количество узлов: " + nodeNum);
        System.out.println("Шаг: " + step);
        for (int i = 1; i <= nodeNum / 2; i++) {
            sum1 += calcFunction(intervalBottom + (2 * i - 1) * step);
        }
        for (int i = 1; i <= nodeNum / 2 - 1; i++) {
            sum2 += calcFunction(intervalBottom + 2 * i * step);
        }
        return (step / 3) * (calcFunction(intervalBottom) + 4 * sum1 + 2 * sum2 + calcFunction(intervalUpper));
    }

    private static int calcNodeNumSimpson(double intervalBottom, double intervalUpper, double accuracy) {
        double fourthDerivativeSupremumArgument = 0.420399;
        double fourthDerivativeSupremum = calcFourthDerivative(fourthDerivativeSupremumArgument);
        Double result = Math.ceil(Math.pow(fourthDerivativeSupremum * Math.pow(intervalUpper - intervalBottom, 5) / (180 * accuracy), 0.25));
        return result.intValue();
    }

    private static double calcFourthDerivative(double x) {
        return 72 * Math.pow(x, 2) * (5 * Math.pow(x, 6) - 17 * Math.pow(x, 3) + 5) / Math.pow(Math.pow(x, 3) + 1, 5);
    }

    private static double calcIntegralGauss(double intervalBottom, double intervalUpper, int nodeNum) {
        double[] nodes = {-0.86113631, -0.33998104, 0.33998104, 0.86113631};
        double[] coefs = {0.34785484, 0.65214516, 0.65214516, 0.34785484};
        double sum = 0;
        double t;
        System.out.println("Количество узлов: " + nodeNum);
        for (int i = 0; i < nodeNum; i++) {
            t = ((intervalUpper - intervalBottom) * nodes[i] + intervalBottom + intervalUpper) / 2;
            sum += coefs[i] * calcFunction(t);
        }
        return (intervalUpper - intervalBottom) / 2 * sum;
    }
}

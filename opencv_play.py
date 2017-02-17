import numpy as np
import cv2


image_name = 'test.png'

# Open image
im = cv2.imread(image_name, cv2.IMREAD_GRAYSCALE)

ret, thresh = cv2.threshold(im, 250, 255, cv2.THRESH_BINARY)
cv2.imshow('contours', thresh)
cv2.waitKey(0)
cv2.destroyAllWindows()

im2, contours, hierarchy = cv2.findContours(thresh, cv2.RETR_LIST, 2)
area = [cv2.contourArea(cnt) for cnt in contours]
print(area)

im_with_contours = im
for i, c in enumerate(contours):
    if cv2.contourArea(c) > 10000:
        print(cv2.contourArea(c))
        im_with_contours = cv2.drawContours(im_with_contours, c, -1, (0,1,132), 1)


cv2.imwrite('contours.png', im_with_contours)
cv2.imshow('contours', im_with_contours)
cv2.waitKey(0)
cv2.destroyAllWindows()

# Blob detector play
'''
# Do params
params = cv2.SimpleBlobDetector_Params()

# Change thresholds
params.minThreshold = 0
params.maxThreshold = 255
params.thresholdStep = 254
# params.blobColor = 255
params.minDistBetweenBlobs = .0001

# Filter by Area.
params.filterByArea = True
params.minArea = 1500
params.maxArea = 1e20

# Filter by Circularity
params.filterByCircularity = False
params.minCircularity = 0.0001

# Filter by Convexity
params.filterByConvexity = False
params.minConvexity = 0.1

# Filter by Inertia
params.filterByInertia = False
params.minInertiaRatio = 0.01

# Create detector
detector = cv2.SimpleBlobDetector_create(params)

# Detect blobs.
keypoints = detector.detect(im)
print(keypoints)

# Draw detected blobs as red circles.
# cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS ensures
# the size of the circle corresponds to the size of blob

im_with_keypoints = cv2.drawKeypoints(im, keypoints, np.array([]), (0, 0, 255), cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)

# Show blobs
cv2.imshow("Keypoints", im_with_keypoints)
cv2.waitKey(0)
cv2.destroyAllWindows()
'''



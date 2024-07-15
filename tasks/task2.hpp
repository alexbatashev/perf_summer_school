#include <array>
#include <cstddef>

class EventID {
public:
  EventID() : mId(-1) {}

private:
  friend class EventDispatcher;

  EventID(unsigned char id) : mId(id) {}

  unsigned char mId;
};

class EventDispatcher {
public:
  EventDispatcher() = default;

  EventID acquire() {
    for (std::size_t i = 0; i < mFreeList.size(); i++) {
      if (mFreeList[i] == false) {
        mFreeList[i] = true;
        return EventID(i);
      }
    }

    return EventID();
  }

  void release(EventID id) {
    mFreeList[id.mId] = false;
  }

  bool hasFreeSlots() {
    for (bool b : mFreeList) {
      if (b == false) return true;
    }

    return false;
  }
private:
  std::array<bool, 64> mFreeList = {false};
};

